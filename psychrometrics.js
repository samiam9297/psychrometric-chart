// Inputs
const forecast = {
    0: {
        'date': 'Today',
        'dryBulb': 95,
        'humidityRatio': 0.008
    },
    1: {
        'date': 'Feb 11',
        'dryBulb': 90,
        'humidityRatio': 0.007
    },
    2: {
        'date': 'Feb 12',
        'dryBulb': 85,
        'humidityRatio': 0.006
    },
    3: {
        'date': 'Feb 13',
        'dryBulb': 80,
        'humidityRatio': 0.005
    },
    4: {
        'date': 'Feb 14',
        'dryBulb': 75,
        'humidityRatio': 0.004
    },
    5: {
        'date': 'Feb 15',
        'dryBulb': 70,
        'humidityRatio': 0.003
    },
    6: {
        'date': 'Feb 16',
        'dryBulb': 65,
        'humidityRatio': 0.002
    }
}

/* global ko, d3 */
/* global Blob */
const c8 = -1.0440397e4
const c9 = -1.129465e1
const c10 = -2.7022355e-2
const c11 = 1.289036e-5
const c12 = -2.4780681e-9
const c13 = 6.5459673
const minTemp = 20
const Rda = 53.35 // Dry air gas constant, ft-lbf / lbda-R

function newtonRaphson (zeroFunc, derivativeFunc, initialX, tolerance) {
    if (typeof tolerance === 'undefined') tolerance = 0.0001

    let testX = initialX
    while (Math.abs(zeroFunc(testX)) > tolerance) {
        testX = testX - zeroFunc(testX) / derivativeFunc(testX)
    }
    return testX
}

// Utility method that guarantees that min and max are exactly
// as input, with the step size based on 0.
function range (min, max, stepsize) {
    const parsedMin = parseFloat(min)
    const toReturn = parsedMin % stepsize === 0 ? [] : [parsedMin]
    let n = 0
    const baseValue = stepsize * Math.ceil(parsedMin / stepsize)
    while (baseValue + n * stepsize < parseFloat(max)) {
        toReturn.push(baseValue + n * stepsize)
        n = n + 1
    }

    toReturn.push(max)
    return toReturn
}

// Saturation pressure in psi from temp in °F.
function satPressFromTempIp (temp) {
    const t = temp + 459.67
    const lnOfSatPress =
        c8 / t +
        c9 +
        c10 * t +
        c11 * Math.pow(t, 2) +
        c12 * Math.pow(t, 3) +
        c13 * Math.log(t)
    return Math.exp(lnOfSatPress)
}

function satHumidRatioFromTempIp (temp, totalPressure) {
    if (arguments.length !== 2) throw Error(`Not all parameters specified. temp: ${temp}; P: ${totalPressure}`)
    const satPress = satPressFromTempIp(temp)
    return (0.621945 * satPress) / (totalPressure - satPress)
}

function wFromPv (pv, totalPressure) {
    if (arguments.length !== 2) throw Error(`Not all parameters specified. pv: ${pv}; P: ${totalPressure}`)
    return (0.621945 * pv) / (totalPressure - pv)
}

function pvFromw (w, totalPressure) {
    if (typeof w === 'string') w = parseFloat(w)
    if (w < 0.000001) return 0
    return totalPressure / (1 + 0.621945 / w)
}

// partial pressure of vapor from dry bulb temp (°F) and rh (0-1)
function pvFromTempRh (temp, rh) {
    if (rh < 0 || rh > 1) throw new Error('RH value must be between 0-1')
    return rh * satPressFromTempIp(temp)
}

function tempFromRhAndPv (rh, pv) {
    if (!rh || rh > 1) throw new Error('RH value must be between 0-1')

    const goalPsat = pv / rh

    // Employ Newton-Raphson method.
    function funcToZero (temp) {
        return satPressFromTempIp(temp) - goalPsat
    }

    const derivativeFunc = (temp) => dPvdT(1, temp)
    return newtonRaphson(funcToZero, derivativeFunc, 80, 0.00001)
}

function tempFromEnthalpyPv (h, pv, totalPressure) {
    const w = wFromPv(pv, totalPressure)
    return (h - w * 1061) / (0.24 + w * 0.445)
}

// Returns object with temperature (°F) and vapor pressure (psia)
function tempPvFromvRh (v, rh, totalPressure) {
    const rAir = 53.35 // Gas constant in units of ft-lbf / lbm - R

    function funcToZero (temp) {
        // The 144 is a conversion factor from psf to psi. The 469.67 is to go from F to R.
        const term1 = satPressFromTempIp(temp) * rh
        const term2 = (totalPressure - rAir * (temp + 459.67) / (v * 144))
        return term1 - term2
    }

    function derivative (temp) {
        return dPvdT(rh, temp) + rAir / (v * 144)
    }

    // Employ the Newton-Raphson method.
    testTemp = newtonRaphson(funcToZero, derivative, 80)
    return { temp: testTemp, pv: pvFromTempRh(testTemp, rh) }
}

function WetBulbRh (wetBulb, rh, totalP) {
    if (rh < 0 || rh > 1) {
        throw new Error('RH expected to be between 0 and 1')
    }

    function funcToZero (testTemp) {
        const w1 = wFromWetbulbDryBulb(wetBulb, testTemp, totalP)
        const pv2 = rh * satPressFromTempIp(testTemp)
        const w2 = wFromPv(pv2, totalP)
        return w1 - w2
    }

    let updatedMaxTemp = 200
    let updatedMinTemp = 0
    let looping = true

    while (looping) {
        const testTemp = (updatedMaxTemp + updatedMinTemp) / 2

        const result = funcToZero(testTemp)

        if (Math.abs(result) < 0.00001) {
            looping = false
        } else {
            // Too low case
            if (result > 0) {
                updatedMinTemp = testTemp
            } else { updatedMaxTemp = testTemp }
        }
    }

    return { temp: testTemp, pv: pvFromTempRh(testTemp, rh) }
}

// temp: Dry bulb temperature in °F
// w: Humidity ratio
// totalPressure: Total Pressure in psia.
function wetBulbFromTempw (temp, w, totalPressure) {
    // Function we'd like to 0. A difference in w's.
    function testWetbulbResult (testWetbulb) {
        const satwAtWetBulb = satHumidRatioFromTempIp(testWetbulb, totalPressure)

        return ((1093 - 0.556 * testWetbulb) * satwAtWetBulb - 0.24 * (temp - testWetbulb)) /
            (1093 + 0.444 * temp - testWetbulb) - w
    }

    let updatedMaxTemp = temp
    let updatedMinTemp = 0

    let testTemp = (updatedMaxTemp + updatedMinTemp) / 2

    let iterations = 0

    let testResult = testWetbulbResult(testTemp)

    while (Math.abs(testResult) > 0.000001) {
        if (iterations > 500) {
            throw new Error('Infinite loop in temp from Rh and Pv.')
        }

        if (testResult > 0) {
            updatedMaxTemp = testTemp
            testTemp = (updatedMaxTemp + updatedMinTemp) / 2
        } else {
            updatedMinTemp = testTemp
            testTemp = (updatedMaxTemp + updatedMinTemp) / 2
        }

        testResult = testWetbulbResult(testTemp)
        iterations++
    }

    return testTemp
}

function tempFromWetbulbw (wetBulb, w, totalPressure) {
    const wsatWetBulb = satHumidRatioFromTempIp(wetBulb, totalPressure)
    return ((1093 - 0.556 * wetBulb) * wsatWetBulb + 0.24 * wetBulb - w * (1093 - wetBulb)) / (0.444 * w + 0.24)
}

function wFromWetbulbDryBulb (wetbulbTemp, temp, totalPressure) {
    const wsatWetBulb = satHumidRatioFromTempIp(wetbulbTemp, totalPressure)
    return ((1093 - 0.556 * wetbulbTemp) * wsatWetBulb - 0.24 * (temp - wetbulbTemp)) / (1093 + 0.444 * temp - wetbulbTemp)
}

function vFromTempw (temp, w, totalPressure) {
    return 0.370486 * (temp + 459.67) * (1 + 1.607858 * w) / totalPressure
}

function tempFromvw (v, w, totalPressure) {
    return (v * totalPressure) / (0.370486 * (1 + 1.607858 * w)) - 459.67
}

function wFromTempv (temp, v, totalPressure) {
    const numerator = ((totalPressure * v) / (0.370486 * (temp + 459.67))) - 1
    return numerator / 1.607858
}

// Calculate derivative of pv vs. T at given RH (0-1) and temp (°F)
function dPvdT (rh, temp) {
    if (rh < 0 || rh > 1) throw Error('rh should be specified 0-1')
    const absTemp = temp + 459.67
    const term1 =
        -c8 / (absTemp * absTemp) +
        c10 +
        2 * c11 * absTemp +
        3 * c12 * absTemp * absTemp +
        c13 / absTemp
    return rh * satPressFromTempIp(temp) * term1
}

const pixelWidth = 1300
const pixelHeight = 700

const xOffsetPercentLeft = 2
const xOffsetPercentRight = 15
const yOffsetPercent = 10

const yCanvasRange = [
    pixelHeight - (yOffsetPercent * pixelHeight) / 100,
    (yOffsetPercent * pixelHeight) / 100
]

const svg = d3.select('svg')

svg.style('width', `${pixelWidth  }px`)
svg.style('height', `${pixelHeight  }px`)


function humidityRatioFromEnthalpyTemp (enthalpy, temp) {
    return (enthalpy - 0.24 * temp) / (1061 + 0.445 * temp)
}

function enthalpyFromTempPv (temp, pv, totalPressure) {
    const w = wFromPv(pv, totalPressure)
    return 0.24 * temp + w * (1061 + 0.445 * temp)
}

function pvFromEnthalpyTemp (enthalpy, temp, totalPressure) {
    return pvFromw(humidityRatioFromEnthalpyTemp(enthalpy, temp), totalPressure)
}

function satTempAtEnthalpy (enthalpy, totalPressure) {
    let currentLowTemp = 0
    let currentHighTemp = 200

    let error = 1
    let testTemp = (currentLowTemp + currentHighTemp) / 2

    let iterations = 0
    do {
        iterations++
        if (iterations > 500) throw Error('Inifite loop in satTempAtEnthalpy')
        testTemp = (currentLowTemp + currentHighTemp) / 2
        const testSatHumidityRatio = satHumidRatioFromTempIp(testTemp, totalPressure)
        const testHumidityRatio = humidityRatioFromEnthalpyTemp(
            enthalpy,
            testTemp
        )

        error = testSatHumidityRatio - testHumidityRatio
        if (testSatHumidityRatio > testHumidityRatio) {
            currentHighTemp = testTemp
        } else {
            currentLowTemp = testTemp
        }
    } while (Math.abs(error) > 0.00005)

    return testTemp
}

function isMult (val, mult) { return val % mult === 0 }

const constantRHvalues = [10, 20, 30, 40, 50, 60, 70, 80, 90]

function StateTempw (maxTemp, maxw, index, totalPressure) {
    const self = this

    self.temperature = ko.observable(forecast[index].dryBulb)
    self.humidityRatio = ko.observable(forecast[index].humidityRatio)

    self.pv = ko.computed(() => pvFromw(self.humidityRatio(), totalPressure))
    self.name = ko.observable(forecast[index].date)
}

function ViewModel () {
    const self = this
    // Start by creating svg elements in the order that I want
    // them layered. The later items will be on top of the earlier items.
    svg.append('g').attr('id', 'modules-1')
    svg.append('g').attr('id', 'modules-2')
    svg.append('g').attr('id', 'modules-3')
    svg.append('g').attr('id', 'modules-4')
    svg.append('g').attr('id', 'modules-1-legend')
    svg.append('g').attr('id', 'modules-2-legend')
    svg.append('g').attr('id', 'modules-3-legend')
    svg.append('g').attr('id', 'modules-4-legend')
    svg.append('g').attr('id', 'modules-5-legend')
    svg.append('g').attr('id', 'specific-humidity-lines')
    svg.append('g').attr('id', 'x-axis')
    const wetBulbPaths = svg.append('g').attr('id', 'wetbulb-lines')
    svg.append('g').attr('id', 'yAxisHumid').style('visibility', 'hidden')
    svg.append('g').attr('id', 'rh-lines')
    svg.append('g').attr('id', 'temp-lines')

    svg.append('g').attr('id', 'h-labels')
    svg.append('g').attr('id', 'boundary-lines').append('path')
        .attr('stroke', '#000000')
        .attr('stroke-width', 2)
        .attr('fill', 'none')

    svg.append('g').attr('id', 'rh-label-background')
    const rhticks = svg
        .append('g')
        .attr('class', 'ticks')
        .attr('id', 'rh-ticks')

    svg.append('g').attr('id', 'v-label-backgrounds')
    svg.append('g').attr('id', 'v-labels')

    svg.append('g').attr('id', 'wetbulb-labels-backgrounds')
    svg.append('g').attr('id', 'wetbulb-labels')

    svg.append('g').attr('id', 'modules-labels')
    svg.append('g').attr('id', 'states')
    svg.append('g').attr('id', 'state-circles')
    svg.append('g').attr('id', 'state-backgrounds')
    svg.append('g').attr('id', 'state-text')
    svg.append('g').attr('id', 'dewpointlabels')

    self.maxTempInput = ko.observable('100').extend({ rateLimit: 500 })
    self.maxTemp = ko.computed(() => {
        const parsedValue = parseInt(self.maxTempInput())
        if (!isNaN(parsedValue) && parsedValue > minTemp && parsedValue < 180) return parsedValue
        return 100
    })

    self.totalPressureInput = ko.observable('14.7').extend({ rateLimit: 500 })
    self.totalPressure = ko.pureComputed(() => {
        const parsedValue = parseFloat(self.totalPressureInput())
        if (!isNaN(parsedValue) && parsedValue > 10 && parsedValue < 20) return parsedValue
        return 14.7
    })

    self.maxwInput = ko.observable('0.03').extend({ rateLimit: 500 })
    self.maxw = ko.pureComputed(() => {
        const parsedValue = parseFloat(self.maxwInput())
        if (!isNaN(parsedValue) && parsedValue > 0 && parsedValue < 0.07) return parsedValue
        return 0.03
    })

    self.xScale = ko.computed(() => {
        return d3.scaleLinear()
            .domain([minTemp, self.maxTemp()])
            .range([
                (xOffsetPercentLeft * pixelWidth) / 100,
                pixelWidth - (xOffsetPercentRight * pixelWidth) / 100
            ])
    })

    self.pixelsPerTemp = ko.pureComputed(() => self.xScale()(1) - self.xScale()(0))
    self.pixelsPerPsia = ko.pureComputed(() => self.yScale()(1) - self.yScale()(0))

    // Return angle in °, given slope in units of psi / °F
    angleFromDerivative = derivative => (Math.atan(derivative * self.pixelsPerPsia() / (self.pixelsPerTemp())
    ) * 180) / Math.PI

    self.maxPv = ko.pureComputed(() => pvFromw(self.maxw(), self.totalPressure()))

    self.yScale = ko.pureComputed(() => {
        return d3
            .scaleLinear()
            .domain([0, self.maxPv()])
            .range(yCanvasRange)
    })

    self.yAxis = ko.pureComputed(() => {
        return d3.axisRight().scale(self.yScale())
    })

    self.saturationLine = ko.pureComputed(() => {
        return d3
            .line()
            .x(d => self.xScale()(d.x))
            .y(d => self.yScale()(Math.min(d.y, self.maxPv())))
    })

    self.tempAtCutoff = ko.pureComputed(() => tempFromRhAndPv(1, self.maxPv()))
    self.upperLeftBorderTemp = ko.pureComputed(() => {
        return self.tempAtCutoff() - 0.05 * (self.maxTemp() - minTemp)
    })

    self.bottomLeftBorderPv = ko.pureComputed(() => {
        return satPressFromTempIp(minTemp) + 0.05 * self.maxPv()
    })

    self.constantTemps = ko.pureComputed(() => range(minTemp, self.maxTemp(), 5))

    self.constantTempLines = ko.computed(() => {
        return self.constantTemps().map(temp => {
            return [{ x: temp, y: 0 }, { x: temp, y: satPressFromTempIp(temp) }]
        })
    })

    ko.computed(function () {
        const selection = d3.select('#temp-lines')
            .selectAll('path')
            .data(self.constantTempLines())

        selection
            .enter()
            .append('path')
            .merge(selection)
            .attr('d', d => self.saturationLine()(d))
            .attr('fill', 'none')
            .attr('stroke', '#000000')
            .attr('stroke-width', d => d[0].x % 10 === 0 ? 1 : 0.5)

        selection.exit().remove()
    })

    self.constantHumidities = ko.computed(() => {
        const humidityStep = 0.002
        const constantHumidities = []
        for (let i = humidityStep; i < wFromPv(self.maxPv(), self.totalPressure()); i = i + humidityStep) {
            constantHumidities.push(i)
        }
        return constantHumidities
    })

    self.constantHumidityLines = ko.computed(() => {
        return self.constantHumidities().map(humidity => {
            const pv = pvFromw(humidity, self.totalPressure())
            return [
                {
                    x: pv < satPressFromTempIp(minTemp) ? minTemp : tempFromRhAndPv(1, pv),
                    y: pv
                },
                { x: self.maxTemp(), y: pv }
            ]
        })
    })

    ko.computed(() => {
        const selection = d3.select('#specific-humidity-lines').selectAll('path').data(self.constantHumidityLines())
        selection.enter()
            .append('path')
            .attr('fill', 'none')
            .attr('stroke', 'blue')
            .attr('stroke-width', 0.5)
            .merge(selection)
            .attr('d', d => self.saturationLine()(d))

        selection.exit().remove()
    })

    self.xAxis = ko.computed(() => {
        return d3
            .axisBottom()
            .scale(self.xScale())
            .tickValues(range(minTemp, self.maxTemp(), 5).filter(temp => temp % 5 === 0))
    })

    ko.computed(() => {
        d3.select('#x-axis').attr('transform', `translate(0,${  self.yScale()(-0.005)  })`)

        const axis = self.xAxis()
        d3.select('#x-axis').call(axis)
    })

    self.yAxisHumid = ko.computed(() => {
        return d3
            .axisRight()
            .scale(self.yScale())
            .tickValues(self.constantHumidities().map(w => pvFromw(w, self.totalPressure())))
            .tickFormat(d => wFromPv(d, self.totalPressure()).toFixed(3))
    })

    ko.computed(() => {
        d3.select('#yAxisHumid')
            .attr('transform', `translate(${  self.xScale()(parseInt(self.maxTemp()) + 0.5)  },0)`)
            .call(self.yAxisHumid())
    })

    // Want the temp diff to be 10% of total width, 9 labels.
    const tempdiff = ko.pureComputed(() => Math.round((self.maxTemp() - minTemp) * 0.15 / 9))
    const starttemp = ko.pureComputed(() => Math.round(minTemp + (self.maxTemp() - minTemp) * 0.6))

    self.constRHLines = ko.computed(() => {
        return constantRHvalues.map((rhValue, i) => {
            const mapFunction = temp => ({
                x: temp,
                y: (satPressFromTempIp(temp) * rhValue) / 100
            })
            let data
            if (pvFromTempRh(self.maxTemp(), rhValue / 100) < self.maxPv()) {
                data = range(minTemp, self.maxTemp(), 0.5).map(mapFunction)
            } else {
                const tempAtBorder = tempFromRhAndPv(rhValue / 100, self.maxPv())
                data = range(minTemp, tempAtBorder, 0.5).map(mapFunction)
            }

            const temp = starttemp() - i * tempdiff()
            const pv = pvFromTempRh(temp, rhValue / 100)

            //// Get derivative in psia/°F
            const derivative = dPvdT(rhValue / 100, temp)
            //// Need to get in same units, pixel/pixel
            const rotationDegrees = angleFromDerivative(derivative)

            return {
                rh: rhValue,
                temp,
                pv,
                data,
                rotationDegrees,
                x: self.xScale()(temp),
                y: self.yScale()(pv)
            }
        })
    })

    ko.computed(() => {
        let selection = d3.select('#rh-lines').selectAll('path').data(self.constRHLines())
        selection
            .enter()
            .append('path')
            .attr('fill', 'none')
            .attr('stroke', 'black')
            .attr('stroke-width', 1)
            .merge(selection)
            .attr('d', d => self.saturationLine()(d.data))

        selection.exit().remove()

        const height = 12
        const labelData = self.constRHLines().filter(d => d.pv < self.maxPv())
        selection = d3.select('#rh-label-background').selectAll('rect').data(labelData)
        selection
            .enter()
            .append('rect')
            .attr('width', 25)
            .attr('height', height)
            .attr('fill', 'white')
            .merge(selection)
            .attr('x', d => self.xScale()(d.temp))
            .attr('y', d => self.yScale()(d.pv))
            .attr('transform', d => `rotate(${d.rotationDegrees}, ${d.x}, ${d.y}) translate(-2 -${height + 2})`)
        selection.exit().remove()

        selection = rhticks.selectAll('text').data(labelData)
        selection.enter()
            .append('text')
            .attr('class', 'rh-ticks')
            .text(d => `${d.rh  }%`)
            .merge(selection)
            .attr('x', d => d.x)
            .attr('y', d => d.y)
            .attr('transform', d => `rotate(${d.rotationDegrees}, ${d.x}, ${d.y}) translate(0 -3)`)
        selection.exit().remove()
    })

    self.minv = ko.computed(() => vFromTempw(minTemp, 0, self.totalPressure()))

    self.maxv = ko.computed(() => vFromTempw(self.maxTemp(), wFromPv(self.maxPv(), self.totalPressure()), self.totalPressure()))
    self.vValues = ko.computed(() => range(Math.ceil(self.minv() / 0.1) * 0.1, Math.floor(self.maxv() / 0.1) * 0.1, 0.1))

    self.vLines = ko.computed(() => {

        const firstVCutoff = vFromTempw(minTemp, satHumidRatioFromTempIp(minTemp, self.totalPressure()), self.totalPressure())
        const secondVCutoff = vFromTempw(self.tempAtCutoff(), wFromPv(self.maxPv(), self.totalPressure()), self.totalPressure())

        return self.vValues().map(v => {
            const mapFunction = temp => {
                return {x: temp, y: pvFromw(wFromTempv(temp, v, self.totalPressure()), self.totalPressure())}
            }
            let lowerTemp
            let upperTemp

            if (v < firstVCutoff) {
                lowerTemp = minTemp
                upperTemp = tempFromvw(v, 0, self.totalPressure())
            } else if (v < secondVCutoff) {
                lowerTemp = tempPvFromvRh(v, 1, self.totalPressure()).temp
                upperTemp = Math.min(tempFromvw(v, 0, self.totalPressure()), self.maxTemp())
            } else {
                lowerTemp = tempFromvw(v, wFromPv(self.maxPv(), self.totalPressure()), self.totalPressure())
                upperTemp = Math.min(tempFromvw(v, 0, self.totalPressure()), self.maxTemp())
            }

            const data = [lowerTemp, upperTemp].map(mapFunction)
            const labelLocation = tempPvFromvRh(v, 0.35, self.totalPressure())

            // 144 to go from psf to psi.
            const derivative = -Rda / v / 144
            const rotationDegrees = angleFromDerivative(derivative)

            return {
                v: Math.round(v * 10) / 10, // properly round to 1 decimal place, because Javascript.
                data,
                labelLocation,
                rotationDegrees,
                x: self.xScale()(labelLocation.temp),
                y: self.yScale()(labelLocation.pv)
            }
        })
    })

    function tempAtStraightEnthalpyLine (enthalpy) {
        const rise = self.maxPv() - self.bottomLeftBorderPv()
        const run = (self.upperLeftBorderTemp()) - minTemp

        function straightLinePv (temp) {
            return self.bottomLeftBorderPv() + (rise / run) * (temp - minTemp)
        }

        function funcToZero (temp) {
            return straightLinePv(temp) - pvFromEnthalpyTemp(enthalpy, temp, self.totalPressure())
        }

        // This comes from maxima, a computer algebra system, see corresponding maxima file.
        function derivative (temp) {
            return (rise / run) - ((1807179 * (12000000 * temp - 50000000 * enthalpy) * self.totalPressure()) /
                Math.pow(1807179 * temp + 50000000 * enthalpy + 32994182250, 2) -
                (12000000 * self.totalPressure()) / (1807179 * temp + 50000000 * enthalpy +
                    32994182250))
        }

        return newtonRaphson(funcToZero, derivative, 80)
    }

    self.minEnthalpy = ko.pureComputed(() => enthalpyFromTempPv(minTemp, 0, self.totalPressure()))
    self.maxEnthalpy = ko.pureComputed(() => {
        return enthalpyFromTempPv(self.maxTemp(), self.maxPv(), self.totalPressure())
    })

    self.constEnthalpyValues = ko.pureComputed(() => {
        return range(Math.ceil(self.minEnthalpy()), Math.floor(self.maxEnthalpy()), 0.2)
    })

    self.enthalpyValueToLine = enthalpyValue => {
        const firstBoundaryEnthalpy = enthalpyFromTempPv(minTemp, satPressFromTempIp(minTemp) + 0.05 * self.maxPv(), self.totalPressure())
        const secondBoundaryEnthalpy = enthalpyFromTempPv(self.upperLeftBorderTemp(), self.maxPv(), self.totalPressure())

        const maxEnthalpyTemp = Math.min(enthalpyValue / 0.24, self.maxTemp())
        const mapFunction = temp => {
            return {x: temp, y: pvFromEnthalpyTemp(enthalpyValue, temp, self.totalPressure())}
        }
        if (enthalpyValue < firstBoundaryEnthalpy) {
            if (enthalpyValue % 5 === 0) {
                return { h: enthalpyValue, coords: range(minTemp, maxEnthalpyTemp, 0.25).map(mapFunction) }
            } else {
                return { h: enthalpyValue, coords: range(minTemp, satTempAtEnthalpy(enthalpyValue, self.totalPressure()), 0.25).map(mapFunction) }
            }
        } else if (enthalpyValue < secondBoundaryEnthalpy) {
            const tempAtBorder = tempAtStraightEnthalpyLine(enthalpyValue)
            return { h: enthalpyValue, coords: range(tempAtBorder, enthalpyValue % 5 === 0 ? maxEnthalpyTemp : satTempAtEnthalpy(enthalpyValue, self.totalPressure()), 0.25).map(mapFunction) }
        } else { // Top section
            return { h: enthalpyValue,
                coords: range(tempFromEnthalpyPv(enthalpyValue, self.maxPv(), self.totalPressure()),
                    isMult(enthalpyValue, 5) ? maxEnthalpyTemp : satTempAtEnthalpy(enthalpyValue, self.totalPressure()), 0.25).map(mapFunction)
            }
        }
    }

    const minWetBulb = wetBulbFromTempw(minTemp, 0, self.totalPressure())
    self.maxWetBulb = ko.computed(() => wetBulbFromTempw(self.maxTemp(), wFromPv(self.maxPv(), self.totalPressure()), self.totalPressure()))
    self.wetBulbBottomRight = ko.computed(() => wetBulbFromTempw(self.maxTemp(), 0, self.totalPressure()))
    self.wetBulbValues = ko.computed(() => range(Math.ceil(minWetBulb), Math.floor(self.maxWetBulb()), 5))

    const wetBulbLabelRh = 0.57 // RH value to put all the wetbulb labels.
    self.wetBulbLines = ko.computed(() => {

        // This is the derivative of Pv vs. temperature for a given
        // constant wet-bulb line.
        derivative = (temp, wetbulb) => {
            const wsatwetbulb = satHumidRatioFromTempIp(wetbulb, self.totalPressure())

            const high = (1093 - 0.556 * wetbulb) * wsatwetbulb - 0.24 * (temp - wetbulb)
            const low = 1093 + 0.444 * temp - wetbulb

            const dHigh = -0.24
            const dLow = 0.444

            const dwdT = ((low * dHigh) - (high * dLow)) / (low * low)

            const w = wFromWetbulbDryBulb(wetbulb, temp, self.totalPressure())

            const dpvdw = (200000 * self.totalPressure()) / (200000 * w + 124389) - (40000000000 * self.totalPressure() * w) / Math.pow(200000 * w + 124389, 2)

            return dpvdw * dwdT
        }

        return self.wetBulbValues().map((wetbulbTemp) => {
            const mapFunction = temp => {
                return {
                    y: pvFromw(wFromWetbulbDryBulb(wetbulbTemp, temp, self.totalPressure()), self.totalPressure()),
                    x: temp
                }
            }

            let lowerTemp
            let upperTemp
            if (wetbulbTemp < minTemp) {
                lowerTemp = minTemp
                upperTemp = tempFromWetbulbw(wetbulbTemp, 0, self.totalPressure())
            } else if (wetbulbTemp < self.wetBulbBottomRight()) {
                lowerTemp = wetbulbTemp
                upperTemp = tempFromWetbulbw(wetbulbTemp, 0, self.totalPressure())
            } else if (wetbulbTemp < self.tempAtCutoff()) {
                lowerTemp = wetbulbTemp
                upperTemp = self.maxTemp()
            } else {
                lowerTemp = tempFromWetbulbw(wetbulbTemp, wFromPv(self.maxPv(), self.totalPressure()), self.totalPressure())
                upperTemp = self.maxTemp()
            }

            const data = range(lowerTemp, upperTemp, 1).map(mapFunction)
            const labelState = WetBulbRh(wetbulbTemp, wetBulbLabelRh, self.totalPressure())
            const midtemp = labelState.temp
            const rotationAngle = angleFromDerivative(derivative(midtemp, wetbulbTemp))
            const midpv = labelState.pv

            return {
                wetbulbTemp,
                data,
                midtemp,
                midpv,
                x: self.xScale()(midtemp),
                y: self.yScale()(midpv),
                rotationAngle
            }
        })
    })

    // Drawing wet-bulb items.
    ko.computed(() => {
        let selection = wetBulbPaths.selectAll('path').data(self.wetBulbLines())
        selection.enter().append('path')
            .attr('fill', 'none')
            .attr('stroke', 'black')
            .attr('stroke-dasharray', '6')
            .attr('stroke-width', 1)
            .merge(selection)
            .attr('d', d => self.saturationLine()(d.data))
        selection.exit().remove()

        const data = self.wetBulbLines().filter(d => d.wetbulbTemp % 5 === 0 && d.midtemp > minTemp && d.midtemp < self.maxTemp() && d.midpv < self.maxPv())
        selection = d3.select('#wetbulb-labels').selectAll('text').data(data)
        selection.enter()
            .append('text')
            .attr('class', 'ticks')
            .style('font-size', '8px')
            .text(d => d.wetbulbTemp.toFixed(0))
            .merge(selection)
            .attr('x', d => self.xScale()(d.midtemp))
            .attr('y', d => self.yScale()(d.midpv))
            .attr('transform', d => `rotate(${d.rotationAngle}, ${d.x}, ${d.y}) translate(0 -3)`)

        selection.exit().remove()

        selection = d3.select('#wetbulb-labels-backgrounds').selectAll('rect').data(data)
        selection.enter()
            .append('rect')
            .attr('fill', 'white')
            .attr('width', '14px')
            .attr('height', '10px')
            .merge(selection)
            .attr('x', d => self.xScale()(d.midtemp))
            .attr('y', d => self.yScale()(d.midpv))
            .attr('transform', d => `rotate(${d.rotationAngle}, ${d.x}, ${d.y}) translate(0 -3) translate(-2 -8)`)
        selection.exit().remove()
    })

    // Modules Required Area Charts
    ko.computed(() => {
        const wetbulbTemps = [55, 63.5, 70.7, 77, 80.1]
        const moduleBoundsData = []
        wetbulbTemps.forEach((wetbulbTemp) => {
            const mapFunction = temp => {
                return {
                    y: pvFromw(wFromWetbulbDryBulb(wetbulbTemp, temp, self.totalPressure()), self.totalPressure()),
                    x: temp
                }
            }

            let lowerTemp
            let upperTemp
            if (wetbulbTemp < minTemp) {
                lowerTemp = minTemp
                upperTemp = tempFromWetbulbw(wetbulbTemp, 0, self.totalPressure())
            } else if (wetbulbTemp < self.wetBulbBottomRight()) {
                lowerTemp = wetbulbTemp
                upperTemp = tempFromWetbulbw(wetbulbTemp, 0, self.totalPressure())
            } else if (wetbulbTemp < self.tempAtCutoff()) {
                lowerTemp = wetbulbTemp
                upperTemp = self.maxTemp()
            } else {
                lowerTemp = tempFromWetbulbw(wetbulbTemp, wFromPv(self.maxPv(), self.totalPressure()), self.totalPressure())
                upperTemp = self.maxTemp()
            }

            const data = range(lowerTemp, upperTemp, 1).map(mapFunction)
            moduleBoundsData.push(data)
        })

        const areaData = [
            [
                {
                    'x0': 55,
                    'y0': 0.21410075364760797,
                    'x1': 63.5,
                    'y1': 0.29011100478567914
                },
                {
                    'x0': 95,
                    'y0': 0.003633448333167196,
                    'x1': 95,
                    'y1': 0.124235434
                },
                {
                    'x0': 95.6920634038409,
                    'y0': 0,
                    'x1': 95.6920634038409,
                    'y1': 0.124235434
                },
                {
                    'x0': 100,
                    'y0': 0,
                    'x1': 100,
                    'y1': 0.02451059
                }
            ],
            [
                {
                    'x0': 63.5,
                    'y0': 0.29011100478567914,
                    'x1': 70.7,
                    'y1': 0.37204468718892597
                },
                {
                    'x0': 95,
                    'y0': 0.1242962844118115,
                    'x1': 95,
                    'y1': 0.24421792322919358
                },
                {
                    'x0': 95.6920634038409,
                    'y0': 0.124235434,
                    'x1': 95.6920634038409,
                    'y1': 0.24421792322919358
                },
                {
                    'x0': 100,
                    'y0': 0.02451059,
                    'x1': 100,
                    'y1': 0.132723886
                }
            ],
            [
                {
                    'x0': 70.7,
                    'y0': 0.37204468718892597,
                    'x1': 77,
                    'y1': 0.4596557939864845
                },
                {
                    'x0': 95,
                    'y0': 0.24421792322919358,
                    'x1': 95,
                    'y1': 0.3651
                },
                {
                    'x0': 95.6920634038409,
                    'y0': 0.24421792322919358,
                    'x1': 95.6920634038409,
                    'y1': 0.3651
                },
                {
                    'x0': 100,
                    'y0': 0.132723886,
                    'x1': 100,
                    'y1': 0.2416
                }
            ],
            [
                {
                    'x0': 77,
                    'y0': 0.4596557939864845,
                    'x1': 80.1,
                    'y1': 0.5090155347971046
                },
                {
                    'x0': 95,
                    'y0': 0.3651,
                    'x1': 95,
                    'y1': 0.430822625
                },
                {
                    'x0': 95.6920634038409,
                    'y0': 0.3651,
                    'x1': 95.6920634038409,
                    'y1': 0.430822625

                },
                {
                    'x0': 100,
                    'y0': 0.2416,
                    'x1': 100,
                    'y1': 0.300648804
                }
            ]
        ]

        const colors = ['#0085ac', '#77777a', '#004069', '#a62337']
        const svgYDiv = [9, 4.5, 2.8, 1.9]
        const text = ['1 Module', '1.5 Modules', '2 Modules', '2.5 Modules']

        areaData.forEach((area, index) => {
            const selection = d3.select(`#modules-${  index + 1}`).selectAll('path').data(area)

            const areaFunc = d3.area().curve(d3.curveBasis)
                .x0(function (d) {
                    return self.xScale()(d.x0)
                })
                .x1(function (d) {
                    return self.xScale()(d.x1)
                })
                .y0(function (d) {
                    return self.yScale()(d.y0)
                })
                .y1(function (d) {
                    return self.yScale()(d.y1)
                })

            selection.enter()
                .append('path')
                .attr('fill', colors[index])
                .merge(selection)
                .attr('d', areaFunc(area))
                .style('opacity', 0.65)
            selection.exit().remove()

            // For some reason, D3JS draws duplicate area paths for each entry and overlays them.
            // This is how we'll clean up the duplicates and only keep the first one..
            const cleanUpDiv = document.getElementById(`modules-${  index + 1}`)
            while (cleanUpDiv.children.length > 1) {
                cleanUpDiv.removeChild(cleanUpDiv.lastChild)
            }

            // Module Labels
            const labelDiv = d3.select('#modules-labels')
            const x = self.xScale()(self.maxTemp() + 3)
            const y = self.yScale()(self.maxPv() / svgYDiv[index])
            const width = 100
            const height = 35

            const label = labelDiv.append('g')
                .attr('x', x)
                .attr('y', y)

            label.append('rect')
                .attr('width', width)
                .attr('height', height)
                .attr('fill', colors[index])
                .attr('x', x)
                .attr('y', y)

            label.append('polygon')
                .attr('fill', colors[index])
                .attr('points', [
                    [self.xScale()(self.maxTemp()) - 20, y + 45],
                    [x + 1, y + height / 2 - 5],
                    [x + 1, y + height / 2 + 5]
                ])

            label.append('text')
                .attr('id', `module-${  index + 1  }-label`)
                .attr('font-family', 'FranklinGothic-MediumCond')
                .attr('fill', 'white')
                .text(text[index])
                .attr('x', x + width / 2)
                .attr('y', y + height / 1.6)
                .attr('text-anchor', 'middle')
        })

        const legendLines = [
            [
                {
                    'x': 55,
                    'y': 0.21410075364760797
                },
                {
                    'x': 95.6920634038409,
                    'y': 0
                }
            ],
            [
                {
                    'x': 63.5,
                    'y': 0.29011100478567914
                },
                {
                    'x': 95,
                    'y': 0.1242962844118115
                },
                {
                    'x': 95.6920634038409,
                    'y': 0.124235434
                },
                {
                    'x': 100,
                    'y': 0.02451059
                }
            ],
            [
                {
                    'x': 70.7,
                    'y': 0.37204468718892597
                },
                {
                    'x': 95,
                    'y': 0.24421792322919358
                },
                {
                    'x': 95.6920634038409,
                    'y': 0.24421792322919358
                },
                {
                    'x': 100,
                    'y': 0.132723886
                }
            ],
            [
                {
                    'x': 77,
                    'y': 0.4596557939864845
                },
                {
                    'x': 95,
                    'y': 0.3651
                },
                {
                    'x': 95.6920634038409,
                    'y': 0.3651
                },
                {
                    'x': 100,
                    'y': 0.2416
                }
            ],
            [
                {
                    'x': 80.1,
                    'y': 0.5090155347971046
                },
                {
                    'x': 95,
                    'y': 0.430822625
                },
                {
                    'x': 95.6920634038409,
                    'y': 0.430822625
                },
                {
                    'x': 100,
                    'y': 0.300648804
                }
            ]
        ]

        legendLines.forEach((legend, index) => {
            const selection = d3.select(`#modules-${  index + 1  }-legend`).selectAll('path').data(legend)

            const lineFunc = d3.line().curve(d3.curveBasis)
                .x(function (d) {
                    return self.xScale()(d.x)
                })
                .y(function (d) {
                    return self.yScale()(d.y)
                })

            selection.enter()
                .append('path')
                .attr('fill', 'none')
                .attr('stroke', 'black')
                .attr('stroke-width', 4)
                .attr('stroke-dasharray', index < 3 ? '7' : '0')
                .attr('d', lineFunc(legend))
            selection.exit().remove()
        })
    })
    // END Modules Required Area

    self.boundaryLineData = ko.computed(() => {
        return [
            { x: self.maxTemp(), y: 0 },
            { x: minTemp, y: 0 },
            { x: minTemp, y: satPressFromTempIp(minTemp) },
            ...range(minTemp, tempFromRhAndPv(1, self.maxPv()), 0.1).map((temp) => { return { x: temp, y: satPressFromTempIp(temp) } }),
            { x: tempFromRhAndPv(1, self.maxPv()), y: self.maxPv() },
            { x: self.maxTemp(), y: satPressFromTempIp(tempFromRhAndPv(1, self.maxPv())) },
            { x: self.maxTemp(), y: 0 }
        ]
    })

    ko.computed(() => {
        d3.select('#boundary-lines').select('path')
            .attr('d', `${self.saturationLine()(self.boundaryLineData())  } Z`)
    })

    const statesArray = []
    Object.keys(forecast).forEach((key, index) => {
        statesArray.push(new StateTempw(self.maxTemp(), self.maxw(), index, self.totalPressure()))
    })

    self.states = ko.observableArray(
        statesArray
    )

    self.removeState = (state) => { self.states.remove(state) }

    const elementObservables = [
        {obs: 'showEnthalpyLines', ids: ['enthalpyLines']},
        {obs: 'showw', ids: ['specific-humidity-lines']},
        {obs: 'showTemp', ids: ['temp-lines']},
        {obs: 'showWetBulb', ids: ['wetbulb-lines', 'wetbulb-labels', 'wetbulb-labels-backgrounds']},
        {obs: 'showRh', ids: ['rh-lines', 'rh-ticks', 'rh-label-background']}
    ]

    elementObservables.map(o => {
        if (o.obs !== 'showw') {
            self[o.obs] = ko.observable(true)
        } else {
            self[o.obs] = ko.observable(false)
        }
        ko.computed(() => {
            for (let i = 0; i < o.ids.length; i++) {
                const element = document.getElementById(o.ids[i])
                if (element) {
                    element.style.visibility = self[o.obs]()
                        ? 'visible'
                        : 'hidden'
                }
            }
        })
    })

    ko.computed(() => {
        const rightOffset = 10

        let selection = d3.select('#state-text').selectAll('text').data(self.states())
        selection
            .enter()
            .append('text')
            .merge(selection)
            .attr('x', d => self.xScale()(d.temperature()))
            .attr('y', d => self.yScale()(d.pv()))
            .attr('dx', rightOffset)
            .attr('dy', '-10')
            .text((d) => d.name())
        selection.exit().remove()

        // Once the text has been created we can get the
        // the size of the bounding box to put the background
        // behind.
        const boundingBoxes = []
        d3.select('#state-text').selectAll('text').each(function (d, i) {
            boundingBoxes[i] = this.getBoundingClientRect()
        })

        selection = d3.select('#state-backgrounds').selectAll('rect').data(self.states())
        selection
            .enter()
            .append('rect')
            .merge(selection)
            .attr('x', d => self.xScale()(d.temperature()))
            .attr('y', d => self.yScale()(d.pv()))
            .attr('transform', (d, i) => `translate(${rightOffset - Math.round(boundingBoxes[i].width * 0.1 / 2)}, -25)`)
            .attr('width', (d, i) => `${Math.ceil(boundingBoxes[i].width * 1.1)}px`)
            .attr('height', '20px')
            .attr('fill', 'white')
        selection.exit().remove()

        selection = d3.select('#state-circles').selectAll('circle').data(self.states())
        selection
            .enter()
            .append('circle')
            .style('fill', '#FF7900')
            .attr('stroke', '#000000')
            .attr('stroke-width', 2)
            .attr('r', '5')
            .merge(selection)
            .attr('cx', d => self.xScale()(d.temperature()))
            .attr('cy', d => self.yScale()(d.pv()))
        selection.exit().remove()
    })

    const middleX = self.xScale()((self.maxTemp() + minTemp) / 2)

    // X-axis label
    svg.append('text')
        .text('DRY BULB TEMPERATURE [°F]')
        .attr('x', middleX)
        .attr('y', self.yScale()(-0.06))
        .attr('font-family', 'FranklinGothic-MediumCond')

    // w label
    svg.append('text')
        .attr('id', 'wlabel')
        .attr('font-family', 'FranklinGothic-MediumCond')
        .text('HUMIDITY RATIO')
        .attr('x', self.xScale()(self.maxTemp() + 4))
        .attr('y', self.yScale()(self.maxPv() / 2))
        .style('visibility', 'hidden')

    // Main enthalpy axis label
    svg.append('text')
        .attr('id', 'enthalpy-label')
        .attr('font-family', 'FranklinGothic-MediumCond')
        .text('WET BULB TEMPERATURE [°F]')

    ko.computed(() => {
        const rise = self.maxPv() - self.bottomLeftBorderPv()
        const run = self.upperLeftBorderTemp() - minTemp

        const angle = Math.atan((rise * self.pixelsPerPsia()) / (run * self.pixelsPerTemp())) * 220 / Math.PI

        const basex = 65
        const basey = 0.3

        const absBasex = self.xScale()(basex)
        const absBasey = self.yScale()(basey)


        d3.select('#enthalpy-label')
            .attr('x', absBasex)
            .attr('y', absBasey)
            .attr('transform', `rotate(${angle}, ${absBasex}, ${absBasey}) translate(-100 -40)`)
    })

    // Main enthalpy axis label
    ko.computed(() => {
        const selection = d3.select('#dewpointlabels').selectAll('text')
            .data(
                self.constantTemps().filter(temp => temp % 5 === 0 && satPressFromTempIp(temp) < self.maxPv())
            )
        selection.enter()
            .append('text')
            .text(d => d.toString())
            .attr('dx', '-13')
            .attr('font-size', '10px')
            .merge(selection)
            .attr('x', d => self.xScale()(d))
            .attr('y', d => self.yScale()(satPressFromTempIp(d) + 0.003))
        selection.exit().remove()
    })

    self.blobUrl = ko.pureComputed(() => {
        const blob = new Blob([d3.select('#vizcontainer').node().innerHTML], {type: 'image/svg+xml'})
        return URL.createObjectURL(blob)
    })

}

const viewModel = new ViewModel()
