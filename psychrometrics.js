/* global ko, d3 */
/* global Blob */
/* global saveSvgAsPng */
const c8 = -1.0440397e4;
const c9 = -1.129465e1;
const c10 = -2.7022355e-2;
const c11 = 1.289036e-5;
const c12 = -2.4780681e-9;
const c13 = 6.5459673;

const minTemp = 20;

const Rda = 53.35; // Dry air gas constant, ft-lbf / lbda-R

function getRandomInt(min, max) {
    min = Math.ceil(min);
    max = Math.floor(max);
    return Math.floor(Math.random() * (max - min)) + min; //The maximum is exclusive and the minimum is inclusive
}

function getRandomArbitrary(min, max) {
    return Math.random() * (max - min) + min;
}

function newtonRaphson(zeroFunc, derivativeFunc, initialX, tolerance) {
    if (typeof tolerance === "undefined") tolerance = 0.0001;

    var testX = initialX;
    while(Math.abs(zeroFunc(testX)) > tolerance) {
        testX = testX - zeroFunc(testX) / derivativeFunc(testX);
    }
    return testX;
}

// Utility method that guarantees that min and max are exactly
// as input, with the step size based on 0.
function range(min, max, stepsize) {
    var parsedMin = parseFloat(min);
    var toReturn = parsedMin % stepsize === 0 ? [] : [parsedMin];
    var n = 0;
    var baseValue = stepsize * Math.ceil(parsedMin / stepsize);
    while (baseValue + n * stepsize < parseFloat(max)) {
        toReturn.push(baseValue + n * stepsize);
        n = n + 1;
    }

    toReturn.push(max);
    return toReturn;
}

// Saturation pressure in psi from temp in °F.
function satPressFromTempIp(temp) {
    var t = temp + 459.67;
    var lnOfSatPress =
        c8 / t +
        c9 +
        c10 * t +
        c11 * Math.pow(t, 2) +
        c12 * Math.pow(t, 3) +
        c13 * Math.log(t);
    var satPress = Math.exp(lnOfSatPress);
    return satPress;
}

function satHumidRatioFromTempIp(temp, totalPressure) {
    if (arguments.length !== 2) throw Error(`Not all parameters specified. temp: ${temp}; P: ${totalPressure}`);
    var satPress = satPressFromTempIp(temp);
    return (0.621945 * satPress) / (totalPressure - satPress);
}

function wFromPv(pv, totalPressure) {
    if (arguments.length !== 2) throw Error(`Not all parameters specified. pv: ${pv}; P: ${totalPressure}`);
    return (0.621945 * pv) / (totalPressure - pv);
}

function pvFromw(w, totalPressure) {
    if (typeof w === "string") w = parseFloat(w);
    if (w < 0.000001) return 0;
    return totalPressure / (1 + 0.621945 / w);
}

// partial pressure of vapor from dry bulb temp (°F) and rh (0-1)
function pvFromTempRh(temp, rh) {
    if (rh < 0 || rh > 1) throw new Error("RH value must be between 0-1");
    return rh * satPressFromTempIp(temp);
}

function tempFromRhAndPv(rh, pv) {
    if (!rh || rh > 1) throw new Error("RH value must be between 0-1");

    var goalPsat = pv / rh;

    // Employ Newton-Raphson method.
    function funcToZero(temp) {
        return satPressFromTempIp(temp) - goalPsat;
    }

    var derivativeFunc = (temp) => dPvdT(1, temp);
    return newtonRaphson(funcToZero, derivativeFunc, 80, 0.00001)
}

function tempFromEnthalpyPv(h, pv, totalPressure) {
    var ω = wFromPv(pv, totalPressure);
    return (h - ω * 1061) / (0.24 + ω * 0.445);
}

// Returns object with temperature (°F) and vapor pressure (psia)
function tempPvFromvRh(v, rh, totalPressure) {
    var rAir = 53.35; // Gas constant in units of ft-lbf / lbm - R

    function funcToZero(temp) {
        // The 144 is a conversion factor from psf to psi. The 469.67 is to go from F to R.
        var term1 = satPressFromTempIp(temp) * rh;
        var term2 = (totalPressure - rAir * (temp + 459.67) / (v * 144));
        return term1 - term2;
    }

    function derivative(temp) {
        return dPvdT(rh, temp) + rAir / (v * 144);
    }

    // Employ the Newton-Raphson method.
    testTemp = newtonRaphson(funcToZero, derivative, 80);
    return { temp: testTemp, pv: pvFromTempRh(testTemp, rh) };
}

function WetBulbRh(wetBulb, rh, totalP) {
    if (rh < 0 || rh > 1) {
        throw new Error("RH expected to be between 0 and 1");
    }

    function funcToZero(testTemp) {
        ω1 = ωFromWetbulbDryBulb(wetBulb, testTemp, totalP);
        pv2 = rh * satPressFromTempIp(testTemp);
        ω2 = wFromPv(pv2, totalP);
        return ω1 - ω2;
    }

    var updatedMaxTemp = 200;
    var updatedMinTemp = 0;
    var looping = true;

    while (looping) {
        var testTemp = (updatedMaxTemp + updatedMinTemp) / 2;

        var result = funcToZero(testTemp);

        if (Math.abs(result) < 0.00001) {
            looping = false;
        }
        else {
            // Too low case
            if (result > 0) {
                updatedMinTemp = testTemp;
            }
            else { updatedMaxTemp = testTemp; }
        }
    }

    return { temp: testTemp, pv: pvFromTempRh(testTemp, rh) }
}

// Create Module Bounds array
function createModuleBoundsArea(boundsBottom, boundsTop) {
    let data = [];
    boundsBottom.forEach((a, index) => {
        let x0 = a.x;
        let y0 = a.y;

        data.push({
            x0: x0,
            y0: y0,
            x1: boundsTop[index] ? boundsTop[index].x : 100,
            y1: boundsTop[index] ? boundsTop[index].y : 0
        });
    })

    return data;
}

// temp: Dry bulb temperature in °F
// ω: Humidity ratio
// totalPressure: Total Pressure in psia.
function wetBulbFromTempω(temp, ω, totalPressure) {
    // Function we'd like to 0. A difference in ω's.
    function testWetbulbResult(testWetbulb) {
        var satωAtWetBulb = satHumidRatioFromTempIp(testWetbulb, totalPressure);

        return ((1093 - 0.556 * testWetbulb) * satωAtWetBulb - 0.24 * (temp - testWetbulb)) /
            (1093 + 0.444 * temp - testWetbulb) - ω;
    }

    var updatedMaxTemp = temp;
    var updatedMinTemp = 0;

    var testTemp = (updatedMaxTemp + updatedMinTemp) / 2;

    var iterations = 0;

    var testResult = testWetbulbResult(testTemp);

    while (Math.abs(testResult) > 0.000001) {
        if (iterations > 500) {
            throw new Error("Infinite loop in temp from Rh and Pv.");
        }

        if (testResult > 0) {
            updatedMaxTemp = testTemp;
            testTemp = (updatedMaxTemp + updatedMinTemp) / 2;
        } else {
            updatedMinTemp = testTemp;
            testTemp = (updatedMaxTemp + updatedMinTemp) / 2;
        }

        testResult = testWetbulbResult(testTemp);
        iterations++;
    }

    return testTemp;
}

function tempFromWetbulbω(wetBulb, ω, totalPressure) {
    var ωsatWetBulb = satHumidRatioFromTempIp(wetBulb, totalPressure);
    return ((1093 - 0.556 * wetBulb) * ωsatWetBulb + 0.24 * wetBulb - ω * (1093 - wetBulb)) / (0.444 * ω + 0.24);
}

function ωFromWetbulbDryBulb(wetbulbTemp, temp, totalPressure) {
    var ωsatWetBulb = satHumidRatioFromTempIp(wetbulbTemp, totalPressure);
    return ((1093 - 0.556 * wetbulbTemp) * ωsatWetBulb - 0.24 * (temp - wetbulbTemp)) / (1093 + 0.444 * temp - wetbulbTemp);
}

function vFromTempω(temp, ω, totalPressure) {
    return 0.370486 * (temp + 459.67) * (1 + 1.607858 * ω) / totalPressure;
}

function tempFromvω(v, ω, totalPressure) {
    return (v * totalPressure) / (0.370486 * (1 + 1.607858 * ω)) - 459.67;
}

function ωFromTempv(temp, v, totalPressure) {
    var numerator = ((totalPressure * v) / (0.370486 * (temp + 459.67))) - 1;
    return numerator / 1.607858;
}

// Calculate derivative of pv vs. T at given RH (0-1) and temp (°F)
function dPvdT(rh, temp) {
    if (rh < 0 || rh > 1) throw Error("rh should be specified 0-1");
    var absTemp = temp + 459.67;
    var term1 =
        -c8 / (absTemp * absTemp) +
        c10 +
        2 * c11 * absTemp +
        3 * c12 * absTemp * absTemp +
        c13 / absTemp;
    return rh * satPressFromTempIp(temp) * term1;
}

var pixelWidth = 1300;
var pixelHeight = 700;

var xOffsetPercentLeft = 2;
var xOffsetPercentRight = 15;
var yOffsetPercent = 10;

var yCanvasRange = [
    pixelHeight - (yOffsetPercent * pixelHeight) / 100,
    (yOffsetPercent * pixelHeight) / 100
];

function boundaryLine(element) {
    return element
        .attr("fill", "none")
        .attr("stroke", "#000000")
        .attr("stroke-width", 2);
}

var svg = d3.select("svg");

svg.style("width", pixelWidth + "px");
svg.style("height", pixelHeight + "px");


function humidityRatioFromEnthalpyTemp(enthalpy, temp) {
    return (enthalpy - 0.24 * temp) / (1061 + 0.445 * temp);
}

function enthalpyFromTempPv(temp, pv, totalPressure) {
    var ω = wFromPv(pv, totalPressure);
    return 0.24 * temp + ω * (1061 + 0.445 * temp);
}

function pvFromEnthalpyTemp(enthalpy, temp, totalPressure) {
    return pvFromw(humidityRatioFromEnthalpyTemp(enthalpy, temp), totalPressure);
}

function satTempAtEnthalpy(enthalpy, totalPressure) {
    var currentLowTemp = 0;
    var currentHighTemp = 200;

    var error = 1;
    var testTemp = (currentLowTemp + currentHighTemp) / 2;

    var iterations = 0;
    do {
        iterations++;
        if (iterations > 500) throw Error("Inifite loop in satTempAtEnthalpy");
        testTemp = (currentLowTemp + currentHighTemp) / 2;
        var testSatHumidityRatio = satHumidRatioFromTempIp(testTemp, totalPressure);
        var testHumidityRatio = humidityRatioFromEnthalpyTemp(
            enthalpy,
            testTemp
        );

        error = testSatHumidityRatio - testHumidityRatio;
        if (testSatHumidityRatio > testHumidityRatio) {
            currentHighTemp = testTemp;
        } else {
            currentLowTemp = testTemp;
        }
    } while (Math.abs(error) > 0.00005);

    return testTemp;
}

function isMult(val, mult) { return val % mult === 0; }

var constantRHvalues = [10, 20, 30, 40, 50, 60, 70, 80, 90];

function StateTempω(maxTemp, maxω, name, totalPressure) {
    var self = this;

    self.temperature = ko.observable(getRandomInt(minTemp, maxTemp));
    var maxωrange = Math.min(satHumidRatioFromTempIp(self.temperature(), totalPressure), maxω);

    self.humidityRatio = ko.observable(Math.round(getRandomArbitrary(0, maxωrange) / 0.001) * 0.001);
    self.pv = ko.computed(() => pvFromw(self.humidityRatio(), totalPressure));
    self.name = ko.observable(name);
}

function ViewModel() {
    var self = this;
    // Start by creating svg elements in the order that I want
    // them layered. The later items will be on top of the earlier items.
    svg.append("g").attr("id", "modules-1");
    svg.append("g").attr("id", "modules-2");
    svg.append("g").attr("id", "modules-3");
    svg.append("g").attr("id", "modules-4");
    svg.append("g").attr("id", "modules-1-legend");
    svg.append("g").attr("id", "modules-2-legend");
    svg.append("g").attr("id", "modules-3-legend");
    svg.append("g").attr("id", "modules-4-legend");
    svg.append("g").attr("id", "modules-5-legend");
    svg.append("g").attr("id", "specific-humidity-lines");
    svg.append("g").attr("id", "x-axis");
    var wetBulbPaths = svg.append("g").attr("id", "wetbulb-lines");
    svg.append("g").attr("id", "yAxisHumid");
    svg.append("g").attr("id", "rh-lines");
    svg.append("g").attr("id", "temp-lines");

    var hLabels = svg.append("g").attr("id", "h-labels");
    svg.append("g").attr("id", "boundary-lines").append("path")
        .attr("stroke", "#000000")
        .attr("stroke-width", 2)
        .attr("fill", "none");

    svg.append("g").attr("id", "rh-label-background");
    var rhticks = svg
        .append("g")
        .attr("class", "ticks")
        .attr("id", "rh-ticks");

    svg.append("g").attr("id", "v-label-backgrounds");
    svg.append("g").attr("id", "v-labels");

    svg.append("g").attr("id", "wetbulb-labels-backgrounds");
    svg.append("g").attr("id", "wetbulb-labels");

    svg.append("g").attr("id", "modules-labels");
    svg.append("g").attr("id", "states");
    svg.append("g").attr("id", "state-circles");
    svg.append("g").attr("id", "state-backgrounds");
    svg.append("g").attr("id", "state-text");
    svg.append("g").attr("id", "dewpointlabels");

    self.maxTempInput = ko.observable("100").extend({ rateLimit: 500 });
    self.maxTemp = ko.computed(() => {
        var parsedValue = parseInt(self.maxTempInput());
        if (!isNaN(parsedValue) && parsedValue > minTemp && parsedValue < 180) return parsedValue;
        return 100;
    });

    self.totalPressureInput = ko.observable("14.7").extend({ rateLimit: 500 });
    self.totalPressure = ko.pureComputed(() => {
        var parsedValue = parseFloat(self.totalPressureInput());
        if (!isNaN(parsedValue) && parsedValue > 10 && parsedValue < 20) return parsedValue;
        return 14.7;
    });

    self.maxωInput = ko.observable("0.03").extend({ rateLimit: 500 });
    self.maxω = ko.pureComputed(() => {
        var parsedValue = parseFloat(self.maxωInput());
        if (!isNaN(parsedValue) && parsedValue > 0 && parsedValue < 0.07) return parsedValue;
        return 0.03;
    });

    self.xScale = ko.computed(() => {
        return d3.scaleLinear()
            .domain([minTemp, self.maxTemp()])
            .range([
                (xOffsetPercentLeft * pixelWidth) / 100,
                pixelWidth - (xOffsetPercentRight * pixelWidth) / 100
            ]);
    });

    self.pixelsPerTemp = ko.pureComputed(() => self.xScale()(1) - self.xScale()(0));
    self.pixelsPerPsia = ko.pureComputed(() => self.yScale()(1) - self.yScale()(0));

    // Return angle in °, given slope in units of psi / °F
    angleFromDerivative = derivative =>
                (Math.atan(derivative * self.pixelsPerPsia() / (self.pixelsPerTemp())
                ) * 180) / Math.PI;

    self.maxPv = ko.pureComputed(() => pvFromw(self.maxω(), self.totalPressure()) );

    self.yScale = ko.pureComputed(() => {
        return d3
            .scaleLinear()
            .domain([0, self.maxPv()])
            .range(yCanvasRange);
    });

    self.yAxis = ko.pureComputed(() => {
        return d3.axisRight().scale(self.yScale());
    });

    self.saturationLine = ko.pureComputed(() => {
        return d3
            .line()
            .x(d => self.xScale()(d.x))
            .y(d => self.yScale()(Math.min(d.y, self.maxPv())));
    });

    self.tempAtCutoff = ko.pureComputed(() => tempFromRhAndPv(1, self.maxPv()));
    self.upperLeftBorderTemp = ko.pureComputed(() => {
        return self.tempAtCutoff() - 0.05 * (self.maxTemp() - minTemp);
    });

    self.bottomLeftBorderPv = ko.pureComputed(() => {
        return satPressFromTempIp(minTemp) + 0.05 * self.maxPv();
    });

    self.constantTemps = ko.pureComputed(() => range(minTemp, self.maxTemp(), 5));

    self.constantTempLines = ko.computed(() => {
        return self.constantTemps().map(temp => {
            return [{ x: temp, y: 0 }, { x: temp, y: satPressFromTempIp(temp) }];
        });
    });

    ko.computed(function () {
        var selection = d3.select("#temp-lines")
            .selectAll("path")
            .data(self.constantTempLines());

        selection
            .enter()
            .append("path")
            .merge(selection)
            .attr("d", d => self.saturationLine()(d))
            .attr("fill", "none")
            .attr("stroke", "#000000")
            .attr("stroke-width", d => d[0].x % 10 === 0 ? 1 : 0.5);

        selection.exit().remove();
    });

    self.constantHumidities = ko.computed(() => {
        var humidityStep = 0.002;
        var constantHumidities = [];
        for (let i = humidityStep; i < wFromPv(self.maxPv(), self.totalPressure()); i = i + humidityStep) {
            constantHumidities.push(i);
        }
        return constantHumidities;
    });

    self.constantHumidityLines = ko.computed(() => {
        return self.constantHumidities().map(humidity => {
            var pv = pvFromw(humidity, self.totalPressure());
            return [
                {
                    x: pv < satPressFromTempIp(minTemp) ? minTemp : tempFromRhAndPv(1, pv),
                    y: pv
                },
                { x: self.maxTemp(), y: pv }
            ];
        });
    });

    ko.computed(() => {
        var selection = d3.select("#specific-humidity-lines").selectAll("path").data(self.constantHumidityLines());
        selection.enter()
            .append("path")
            .attr("fill", "none")
            .attr("stroke", "blue")
            .attr("stroke-width", 0.5)
            .merge(selection)
            .attr("d", d => self.saturationLine()(d));

        selection.exit().remove();
    });

    self.xAxis = ko.computed(() => {
        return d3
            .axisBottom()
            .scale(self.xScale())
            .tickValues(range(minTemp, self.maxTemp(), 5).filter(temp => temp % 5 === 0));
    });

    ko.computed(() => {
        d3.select("#x-axis").attr("transform", "translate(0," + self.yScale()(-0.005) + ")");

        var axis = self.xAxis();
        d3.select("#x-axis").call(axis);
    });

    self.yAxisHumid = ko.computed(() => {
        return d3
            .axisRight()
            .scale(self.yScale())
            .tickValues(self.constantHumidities().map(ω => pvFromw(ω, self.totalPressure())))
            .tickFormat(d => wFromPv(d, self.totalPressure()).toFixed(3));
    });

    ko.computed(() => {
        d3.select("#yAxisHumid")
            .attr("transform", "translate(" + self.xScale()(parseInt(self.maxTemp()) + 0.5) + ",0)")
            .call(self.yAxisHumid());
    });

    // Want the temp diff to be 10% of total width, 9 labels.
    var tempdiff = ko.pureComputed(() => Math.round((self.maxTemp() - minTemp) * 0.15 / 9));
    var starttemp = ko.pureComputed(() => Math.round(minTemp + (self.maxTemp() - minTemp) * 0.6));

    self.constRHLines = ko.computed(() => {
        return constantRHvalues.map((rhValue, i) => {
            const mapFunction = temp => ({
                x: temp,
                y: (satPressFromTempIp(temp) * rhValue) / 100
            });
            var data;
            if (pvFromTempRh(self.maxTemp(), rhValue / 100) < self.maxPv()) {
                data = range(minTemp, self.maxTemp(), 0.5).map(mapFunction);
            } else {
                var tempAtBorder = tempFromRhAndPv(rhValue / 100, self.maxPv());
                data = range(minTemp, tempAtBorder, 0.5).map(mapFunction);
            }

            var temp = starttemp() - i * tempdiff();
            var pv = pvFromTempRh(temp, rhValue / 100);

            //// Get derivative in psia/°F
            var derivative = dPvdT(rhValue / 100, temp);
            //// Need to get in same units, pixel/pixel
            var rotationDegrees = angleFromDerivative(derivative);

            return {
                rh: rhValue,
                temp: temp,
                pv: pv,
                data: data,
                rotationDegrees: rotationDegrees,
                x: self.xScale()(temp),
                y: self.yScale()(pv)
            };
        });
    });

    ko.computed(() => {
        var selection = d3.select("#rh-lines").selectAll("path").data(self.constRHLines());
        selection
            .enter()
            .append("path")
            .attr("fill", "none")
            .attr("stroke", "black")
            .attr("stroke-width", 1)
            .merge(selection)
            .attr("d", d => self.saturationLine()(d.data));

        selection.exit().remove();

        var height = 12;
        var labelData = self.constRHLines().filter(d => d.pv < self.maxPv());
        selection = d3.select("#rh-label-background").selectAll("rect").data(labelData);
        selection
            .enter()
            .append("rect")
            .attr("width", 25)
            .attr("height", height)
            .attr("fill", "white")
            .merge(selection)
            .attr("x", d => self.xScale()(d.temp))
            .attr("y", d => self.yScale()(d.pv))
            .attr("transform", d => `rotate(${d.rotationDegrees}, ${d.x}, ${d.y}) translate(-2 -${height + 2})`);
        selection.exit().remove();

        selection = rhticks.selectAll("text").data(labelData);
        selection.enter()
            .append("text")
            .attr("class", "rh-ticks")
            .text(d => d.rh + "%")
            .merge(selection)
            .attr("x", d => d.x)
            .attr("y", d => d.y)
            .attr("transform", d => `rotate(${d.rotationDegrees}, ${d.x}, ${d.y}) translate(0 -3)`);
        selection.exit().remove();
    });

    self.minv = ko.computed(() => vFromTempω(minTemp, 0, self.totalPressure()));

    self.maxv = ko.computed(() => vFromTempω(self.maxTemp(), wFromPv(self.maxPv(), self.totalPressure()), self.totalPressure()));
    self.vValues = ko.computed(() => range(Math.ceil(self.minv() / 0.1) * 0.1, Math.floor(self.maxv() / 0.1) * 0.1, 0.1));

    self.vLines = ko.computed(() => {

        var firstVCutoff = vFromTempω(minTemp, satHumidRatioFromTempIp(minTemp, self.totalPressure()), self.totalPressure());
        var secondVCutoff = vFromTempω(self.tempAtCutoff(), wFromPv(self.maxPv(), self.totalPressure()), self.totalPressure());

        return self.vValues().map(v => {
            var mapFunction = temp => { return { x: temp, y: pvFromw(ωFromTempv(temp, v, self.totalPressure()), self.totalPressure()) }; };
            var lowerTemp;
            var upperTemp;

            if (v < firstVCutoff) {
                lowerTemp = minTemp;
                upperTemp = tempFromvω(v, 0, self.totalPressure());
            } else if (v < secondVCutoff) {
                lowerTemp = tempPvFromvRh(v, 1, self.totalPressure()).temp;
                upperTemp = Math.min(tempFromvω(v, 0, self.totalPressure()), self.maxTemp());
            } else {
                lowerTemp = tempFromvω(v, wFromPv(self.maxPv(), self.totalPressure()), self.totalPressure());
                upperTemp = Math.min(tempFromvω(v, 0, self.totalPressure()), self.maxTemp());
            }

            var data = [lowerTemp, upperTemp].map(mapFunction);
            var labelLocation = tempPvFromvRh(v, 0.35, self.totalPressure());

            // 144 to go from psf to psi.
            var derivative = -Rda / v / 144;
            var rotationDegrees = angleFromDerivative(derivative);

            return {
                v: Math.round(v * 10) / 10, // properly round to 1 decimal place, because Javascript.
                data: data,
                labelLocation: labelLocation,
                rotationDegrees: rotationDegrees,
                x: self.xScale()(labelLocation.temp),
                y: self.yScale()(labelLocation.pv)
            };
        });
    });

    function tempAtStraightEnthalpyLine(enthalpy) {
        var rise = self.maxPv() - self.bottomLeftBorderPv();
        var run = (self.upperLeftBorderTemp()) - minTemp;

        function straightLinePv(temp) {
            return self.bottomLeftBorderPv() + (rise / run) * (temp - minTemp);
        }

        function funcToZero(temp) {
            return straightLinePv(temp) - pvFromEnthalpyTemp(enthalpy, temp, self.totalPressure());
        }

        // This comes from maxima, a computer algebra system, see corresponding maxima file.
        function derivative(temp) {
            return (rise / run) - ((1807179 * (12000000 * temp - 50000000 * enthalpy) * self.totalPressure()) /
              Math.pow(1807179 * temp + 50000000 * enthalpy + 32994182250, 2) -
              (12000000 * self.totalPressure()) / (1807179 * temp + 50000000 * enthalpy +
                                                     32994182250));
        }

        return newtonRaphson(funcToZero, derivative, 80);
    }

    self.minEnthalpy = ko.pureComputed(() => enthalpyFromTempPv(minTemp, 0, self.totalPressure()));
    self.maxEnthalpy = ko.pureComputed(() => {
        return enthalpyFromTempPv(self.maxTemp(), self.maxPv(), self.totalPressure());
    });

    self.constEnthalpyValues = ko.pureComputed(() => {
        return range(Math.ceil(self.minEnthalpy()), Math.floor(self.maxEnthalpy()), 0.2);
    });

    self.enthalpyValueToLine = enthalpyValue => {
        var firstBoundaryEnthalpy = enthalpyFromTempPv(minTemp, satPressFromTempIp(minTemp) + 0.05 * self.maxPv(), self.totalPressure());
        var secondBoundaryEnthalpy = enthalpyFromTempPv(self.upperLeftBorderTemp(), self.maxPv(), self.totalPressure());

        var maxEnthalpyTemp = Math.min(enthalpyValue / 0.24, self.maxTemp());
        var mapFunction = temp => { return { x: temp, y: pvFromEnthalpyTemp(enthalpyValue, temp, self.totalPressure()) }; };
        if (enthalpyValue < firstBoundaryEnthalpy) {
            if (enthalpyValue % 5 === 0) {
                return { h: enthalpyValue, coords: range(minTemp, maxEnthalpyTemp, 0.25).map(mapFunction) };
            } else {
                return { h: enthalpyValue, coords: range(minTemp, satTempAtEnthalpy(enthalpyValue, self.totalPressure()), 0.25).map(mapFunction) };
            }
        } else if (enthalpyValue < secondBoundaryEnthalpy) {
            var tempAtBorder = tempAtStraightEnthalpyLine(enthalpyValue);
            return { h: enthalpyValue, coords: range(tempAtBorder, enthalpyValue % 5 === 0 ? maxEnthalpyTemp : satTempAtEnthalpy(enthalpyValue, self.totalPressure()), 0.25).map(mapFunction) };
        } else { // Top section
            return { h: enthalpyValue,
                coords: range(tempFromEnthalpyPv(enthalpyValue, self.maxPv(), self.totalPressure()),
                    isMult(enthalpyValue, 5) ? maxEnthalpyTemp : satTempAtEnthalpy(enthalpyValue, self.totalPressure()), 0.25).map(mapFunction)
            };
        }
    };

    var minWetBulb = wetBulbFromTempω(minTemp, 0, self.totalPressure());
    self.maxWetBulb = ko.computed(() => wetBulbFromTempω(self.maxTemp(), wFromPv(self.maxPv(), self.totalPressure()), self.totalPressure()));
    self.wetBulbBottomRight = ko.computed(() => wetBulbFromTempω(self.maxTemp(), 0, self.totalPressure()));
    self.wetBulbValues = ko.computed(() => range(Math.ceil(minWetBulb), Math.floor(self.maxWetBulb()), 5));

    var wetBulbLabelRh = 0.57; // RH value to put all the wetbulb labels.
    self.wetBulbLines = ko.computed(() => {

        // This is the derivative of Pv vs. temperature for a given
        // constant wet-bulb line.
        derivative = (temp, wetbulb) => {
            var wsatwetbulb = satHumidRatioFromTempIp(wetbulb, self.totalPressure())

            var high = (1093 - 0.556*wetbulb) * wsatwetbulb - 0.24 * (temp - wetbulb);
            var low = 1093 + 0.444 * temp - wetbulb;

            var dHigh = -0.24;
            var dLow = 0.444;

            var dwdT = ((low * dHigh) - (high * dLow)) / (low * low);

            var w = ωFromWetbulbDryBulb(wetbulb, temp, self.totalPressure());

            var dpvdw = (200000*self.totalPressure())/(200000*w+124389)-(40000000000*self.totalPressure()*w)/Math.pow(200000*w+124389,2);

            return dpvdw * dwdT;
        }

        return self.wetBulbValues().map((wetbulbTemp) => {
            var mapFunction = temp => {
                return {
                    y: pvFromw(ωFromWetbulbDryBulb(wetbulbTemp, temp, self.totalPressure()), self.totalPressure()),
                    x: temp
                };
            };

            var lowerTemp;
            var upperTemp;
            if (wetbulbTemp < minTemp) {
                lowerTemp = minTemp;
                upperTemp = tempFromWetbulbω(wetbulbTemp, 0, self.totalPressure());
            } else if (wetbulbTemp < self.wetBulbBottomRight()) {
                lowerTemp = wetbulbTemp;
                upperTemp = tempFromWetbulbω(wetbulbTemp, 0, self.totalPressure());
            } else if (wetbulbTemp < self.tempAtCutoff()) {
                lowerTemp = wetbulbTemp;
                upperTemp = self.maxTemp();
            } else {
                lowerTemp = tempFromWetbulbω(wetbulbTemp, wFromPv(self.maxPv(), self.totalPressure()), self.totalPressure());
                upperTemp = self.maxTemp();
            }

            var data = range(lowerTemp, upperTemp, 1).map(mapFunction);
            var labelState = WetBulbRh(wetbulbTemp, wetBulbLabelRh, self.totalPressure());
            var midtemp = labelState.temp;
            var rotationAngle = angleFromDerivative(derivative(midtemp, wetbulbTemp));
            var midpv = labelState.pv;

            return {
                wetbulbTemp: wetbulbTemp,
                data: data,
                midtemp: midtemp,
                midpv: midpv,
                x: self.xScale()(midtemp),
                y: self.yScale()(midpv),
                rotationAngle: rotationAngle
            };
        });
    });

    // Drawing wet-bulb items.
    ko.computed(() => {
        var selection = wetBulbPaths.selectAll("path").data(self.wetBulbLines());
        selection.enter().append("path")
            .attr("fill", "none")
            .attr("stroke", "black")
            .attr("stroke-dasharray", "6")
            .attr("stroke-width", 1)
            .merge(selection)
            .attr("d", d => self.saturationLine()(d.data));
        selection.exit().remove();

        var data = self.wetBulbLines().filter(d => d.wetbulbTemp % 5 === 0 && d.midtemp > minTemp && d.midtemp < self.maxTemp() && d.midpv < self.maxPv());
        selection = d3.select("#wetbulb-labels").selectAll("text").data(data);
        selection.enter()
            .append("text")
            .attr("class", "ticks")
            .style("font-size", "8px")
            .text(d => d.wetbulbTemp.toFixed(0))
            .merge(selection)
            .attr("x", d => self.xScale()(d.midtemp))
            .attr("y", d => self.yScale()(d.midpv))
            .attr("transform", d => `rotate(${d.rotationAngle}, ${d.x}, ${d.y}) translate(0 -3)`);

        selection.exit().remove();

        selection = d3.select("#wetbulb-labels-backgrounds").selectAll("rect").data(data);
        selection.enter()
            .append("rect")
            .attr("fill", "white")
            .attr("width", "14px")
            .attr("height", "10px")
            .merge(selection)
            .attr("x", d => self.xScale()(d.midtemp))
            .attr("y", d => self.yScale()(d.midpv))
            .attr("transform", d => `rotate(${d.rotationAngle}, ${d.x}, ${d.y}) translate(0 -3) translate(-2 -8)`);
        selection.exit().remove();
    });

    // Modules Required Area Charts
    ko.computed(() => {
        var wetbulbTemps = [55, 63.5, 70.7, 77, 80.1];
        var moduleBoundsData = [];
        wetbulbTemps.forEach((wetbulbTemp) => {
            var mapFunction = temp => {
                return {
                    y: pvFromw(ωFromWetbulbDryBulb(wetbulbTemp, temp, self.totalPressure()), self.totalPressure()),
                    x: temp
                };
            };

            var lowerTemp;
            var upperTemp;
            if (wetbulbTemp < minTemp) {
                lowerTemp = minTemp;
                upperTemp = tempFromWetbulbω(wetbulbTemp, 0, self.totalPressure());
            } else if (wetbulbTemp < self.wetBulbBottomRight()) {
                lowerTemp = wetbulbTemp;
                upperTemp = tempFromWetbulbω(wetbulbTemp, 0, self.totalPressure());
            } else if (wetbulbTemp < self.tempAtCutoff()) {
                lowerTemp = wetbulbTemp;
                upperTemp = self.maxTemp();
            } else {
                lowerTemp = tempFromWetbulbω(wetbulbTemp, wFromPv(self.maxPv(), self.totalPressure()), self.totalPressure());
                upperTemp = self.maxTemp();
            }

            var data = range(lowerTemp, upperTemp, 1).map(mapFunction);
            moduleBoundsData.push(data)
        })

        var areaData = [
            [
                {
                    "x0": 55,
                    "y0": 0.21410075364760797,
                    "x1": 63.5,
                    "y1": 0.29011100478567914
                },
                {
                    "x0": 95,
                    "y0": 0.003633448333167196,
                    "x1": 95,
                    "y1": 0.124235434
                },
                {
                    "x0": 95.6920634038409,
                    "y0": 0,
                    "x1": 95.6920634038409,
                    "y1": 0.124235434
                },
                {
                    "x0": 100,
                    "y0": 0,
                    "x1": 100,
                    "y1": 0.02451059
                }
            ],
            [
                {
                    "x0": 63.5,
                    "y0": 0.29011100478567914,
                    "x1": 70.7,
                    "y1": 0.37204468718892597
                },
                {
                    "x0": 95,
                    "y0": 0.1242962844118115,
                    "x1": 95,
                    "y1": 0.24421792322919358
                },
                {
                    "x0": 95.6920634038409,
                    "y0": 0.124235434,
                    "x1": 95.6920634038409,
                    "y1": 0.24421792322919358
                },
                {
                    "x0": 100,
                    "y0": 0.02451059,
                    "x1": 100,
                    "y1": 0.132723886
                }
            ],
            [
                {
                    "x0": 70.7,
                    "y0": 0.37204468718892597,
                    "x1": 77,
                    "y1": 0.4596557939864845
                },
                {
                    "x0": 95,
                    "y0": 0.24421792322919358,
                    "x1": 95,
                    "y1": 0.3651
                },
                {
                    "x0": 95.6920634038409,
                    "y0": 0.24421792322919358,
                    "x1": 95.6920634038409,
                    "y1": 0.3651
                },
                {
                    "x0": 100,
                    "y0": 0.132723886,
                    "x1": 100,
                    "y1": 0.2416
                }
            ],
            [
                {
                    "x0": 77,
                    "y0": 0.4596557939864845,
                    "x1": 80.1,
                    "y1": 0.5090155347971046
                },
                {
                    "x0": 95,
                    "y0": 0.3651,
                    "x1": 95,
                    "y1": 0.430822625
                },
                {
                    "x0": 95.6920634038409,
                    "y0": 0.3651,
                    "x1": 95.6920634038409,
                    "y1": 0.430822625

                },
                {
                    "x0": 100,
                    "y0": 0.2416,
                    "x1": 100,
                    "y1": 0.300648804
                }
            ]
        ];

        var colors = ['#0085ac', '#77777a', '#004069', '#a62337']
        var svgYDiv = [9, 4.5, 2.8, 1.9];
        var text = ['1 Module', '1.5 Modules', '2 Modules', '2.5 Modules'];

        areaData.forEach((area, index) => {
            var selection = d3.select("#modules-" + (index + 1)).selectAll("path").data(area);

            var areaFunc = d3.area().curve(d3.curveBasis)
                .x0(function (d) { return self.xScale()(d.x0); })
                .x1(function (d) { return self.xScale()(d.x1); })
                .y0(function (d) { return self.yScale()(d.y0); })
                .y1(function (d) { return self.yScale()(d.y1); });

            selection.enter()
                .append("path")
                .attr("fill", colors[index])
                .merge(selection)
                .attr("d", areaFunc(area))
                .style("opacity", 0.65);
            selection.exit().remove();

            // For some reason, D3JS draws duplicate area paths for each entry and overlays them.
            // This is how we'll clean up the duplicates and only keep the first one..
            var cleanUpDiv = document.getElementById("modules-" + (index + 1));
            while (cleanUpDiv.children.length > 1) {
               cleanUpDiv.removeChild(cleanUpDiv.lastChild);
            }

            // Module Labels
            const labelDiv = d3.select("#modules-labels");
            const x = self.xScale()(self.maxTemp() + 3);
            const y = self.yScale()(self.maxPv() / svgYDiv[index]);
            const width = 100;
            const height = 35;

            var label = labelDiv.append("g")
                .attr("x", x)
                .attr("y", y);

            label.append("rect")
                .attr("width", width)
                .attr("height", height)
                .attr("fill", colors[index])
                .attr("x", x)
                .attr("y", y);

            label.append("polygon")
                .attr("fill", colors[index])
                .attr("points", [
                    [self.xScale()(self.maxTemp()) - 20, y + 45],
                    [x + 1, y + height / 2 - 5],
                    [x + 1, y + height / 2 + 5]
                ]);

            label.append("text")
                .attr("id", "module-" + (index + 1) + "-label")
                .attr("font-family", "FranklinGothic-MediumCond")
                .attr("fill", "white")
                .text(text[index])
                .attr("x", x + width / 2)
                .attr("y", y + height / 1.6)
                .attr("text-anchor", "middle");
        })

        var legendLines = [
            [
                {
                    "x": 55,
                    "y": 0.21410075364760797,
                },
                {
                    "x": 95.6920634038409,
                    "y": 0,
                }
            ],
            [
                {
                    "x": 63.5,
                    "y": 0.29011100478567914,
                },
                {
                    "x": 95,
                    "y": 0.1242962844118115,
                },
                {
                    "x": 95.6920634038409,
                    "y": 0.124235434,
                },
                {
                    "x": 100,
                    "y": 0.02451059,
                }
            ],
            [
                {
                    "x": 70.7,
                    "y": 0.37204468718892597,
                },
                {
                    "x": 95,
                    "y": 0.24421792322919358,
                },
                {
                    "x": 95.6920634038409,
                    "y": 0.24421792322919358,
                },
                {
                    "x": 100,
                    "y": 0.132723886,
                }
            ],
            [
                {
                    "x": 77,
                    "y": 0.4596557939864845,
                },
                {
                    "x": 95,
                    "y": 0.3651,
                },
                {
                    "x": 95.6920634038409,
                    "y": 0.3651,
                },
                {
                    "x": 100,
                    "y": 0.2416,
                }
            ],
            [
                {
                    "x": 80.1,
                    "y": 0.5090155347971046,
                },
                {
                    "x": 95,
                    "y": 0.430822625,
                },
                {
                    "x": 95.6920634038409,
                    "y": 0.430822625,
                },
                {
                    "x": 100,
                    "y": 0.300648804,
                }
            ]
        ];

        legendLines.forEach((legend, index) => {
            var selection = d3.select("#modules-" + (index + 1) + "-legend").selectAll("path").data(legend);
            console.log(selection)

            var lineFunc = d3.line().curve(d3.curveBasis)
                .x(function (d) { return self.xScale()(d.x); })
                .y(function (d) { return self.yScale()(d.y); })

            selection.enter()
                .append("path")
                .attr("fill", "none")
                .attr("stroke", "black")
                .attr("stroke-width", 4)
                .attr("stroke-dasharray", index < 3 ? "7" : "0")
                .attr("d", lineFunc(legend))
            selection.exit().remove();
        })
    });
    // END Modules Required Area

    self.boundaryLineData = ko.computed(() => {
        return [
            { x: self.maxTemp(), y: 0 },
            { x: minTemp, y: 0 },
            { x: minTemp, y: satPressFromTempIp(minTemp) },
            ...range(minTemp, tempFromRhAndPv(1, self.maxPv()), 0.1).map((temp) => { return { x: temp, y: satPressFromTempIp(temp) }; }),
            { x: tempFromRhAndPv(1, self.maxPv()), y: self.maxPv() },
            { x: self.maxTemp(), y: satPressFromTempIp(tempFromRhAndPv(1, self.maxPv())) },
            { x: self.maxTemp(), y: 0 }
        ];
    });

    ko.computed(() => {
        d3.select("#boundary-lines").select("path")
            .attr("d", self.saturationLine()(self.boundaryLineData()) + " Z");
    });

    self.states = ko.observableArray(
        [
            new StateTempω(self.maxTemp(), self.maxω(), "State 1", self.totalPressure())
        ]
    );

    self.addState = () => {
        self.states.push(
            new StateTempω(self.maxTemp(), self.maxω(), "State " + (self.states().length + 1), self.totalPressure())
        );
    };

    self.removeState = (state) => { self.states.remove(state); };

    var elementObservables = [
        { obs: "showEnthalpyLines", ids: ["enthalpyLines"] },
        { obs: "showω", ids: ["specific-humidity-lines"] },
        { obs: "showTemp", ids: ["temp-lines"] },
        { obs: "showWetBulb", ids: ["wetbulb-lines", "wetbulb-labels",  "wetbulb-labels-backgrounds"] },
        { obs: "showRh", ids: ["rh-lines", "rh-ticks", "rh-label-background"] }
    ];

    elementObservables.map(o => {
        if (o.obs !== "showω") {
            self[o.obs] = ko.observable(true);
        } else {
            self[o.obs] = ko.observable(false);
        }
        ko.computed(() => {
            for (let i = 0; i < o.ids.length; i++) {
                var element = document.getElementById(o.ids[i]);
                if (element) {
                    element.style.visibility = self[o.obs]()
                        ? "visible"
                        : "hidden";
                }
            }
        });
    });

    ko.computed(() => {
        var rightOffset = 10;

        var selection = d3.select("#state-text").selectAll("text").data(self.states());
        selection
            .enter()
            .append("text")
            .merge(selection)
            .attr("x", d => self.xScale()(d.temperature()))
            .attr("y", d => self.yScale()(d.pv()))
            .attr("dx", rightOffset)
            .attr("dy", "-10")
            .text((d, i) => d.name());
        selection.exit().remove();

        // Once the text has been created we can get the
        // the size of the bounding box to put the background
        // behind.
        var boundingBoxes = [];
        d3.select("#state-text").selectAll("text").each(function (d, i) {
            boundingBoxes[i] = this.getBoundingClientRect();
        });

        selection = d3.select("#state-backgrounds").selectAll("rect").data(self.states());
        selection
            .enter()
            .append("rect")
            .merge(selection)
            .attr("x", d => self.xScale()(d.temperature()))
            .attr("y", d => self.yScale()(d.pv()))
            .attr("transform", (d, i) => `translate(${rightOffset - Math.round(boundingBoxes[i].width * 0.1 / 2)}, -25)`)
            .attr("width", (d, i) => `${Math.ceil(boundingBoxes[i].width * 1.1)}px`)
            .attr("height", "20px")
            .attr("fill", "white");
        selection.exit().remove();

        selection = d3.select("#state-circles").selectAll("circle").data(self.states());
        selection
            .enter()
            .append("circle")
            .style("fill", "#FF7900")
            .attr("stroke", "#000000")
            .attr("stroke-width", 2)
            .attr("r", "5")
            .merge(selection)
            .attr("cx", d => self.xScale()(d.temperature()))
            .attr("cy", d => self.yScale()(d.pv()));
        selection.exit().remove();
    });

    var middleX = self.xScale()((self.maxTemp() + minTemp) / 2);

    // X-axis label
    svg.append("text")
        .text("DRY BULB TEMPERATURE [°F]")
        .attr("x", middleX)
        .attr("y", self.yScale()(-0.06))
        .attr("font-family", "FranklinGothic-MediumCond");

    // ω label
    svg.append("text")
        .attr("id", "ωlabel")
        .attr("font-family", "FranklinGothic-MediumCond")
        .text("HUMIDITY RATIO")
        .attr("x", self.xScale()(self.maxTemp() + 4))
        .attr("y", self.yScale()(self.maxPv() / 2));

    // Main enthalpy axis label
    svg.append("text")
        .attr("id", "enthalpy-label")
        .attr("font-family", "FranklinGothic-MediumCond")
        .text("WET BULB TEMPERATURE [°F]");

    ko.computed(() => {
        var rise = self.maxPv() - self.bottomLeftBorderPv();
        var run = self.upperLeftBorderTemp() - minTemp;

        var angle = Math.atan((rise * self.pixelsPerPsia()) / (run * self.pixelsPerTemp())) * 220 / Math.PI;

        var basex = 65;
        var basey = 0.3;

        var absBasex = self.xScale()(basex)
        var absBasey = self.yScale()(basey)


        d3.select("#enthalpy-label")
            .attr("x", absBasex)
            .attr("y", absBasey)
            .attr("transform", `rotate(${angle}, ${absBasex}, ${absBasey}) translate(-100 -40)`);
    });

    // Main enthalpy axis label
    ko.computed(() => {
        var selection = d3.select("#dewpointlabels").selectAll("text")
            .data(
                self.constantTemps().filter(temp => temp % 5 === 0 && satPressFromTempIp(temp) < self.maxPv())
            );
        selection.enter()
            .append("text")
            .text(d => d.toString())
            .attr("dx", "-13")
            .attr("font-size", "10px")
            .merge(selection)
            .attr("x", d => self.xScale()(d))
            .attr("y", d => self.yScale()(satPressFromTempIp(d) + 0.003));
        selection.exit().remove();
    });

    self.blobUrl = ko.pureComputed(() => {
        var blob = new Blob([d3.select("#vizcontainer").node().innerHTML], { type: "image/svg+xml" });
        return URL.createObjectURL(blob);
    });

}

var viewModel = new ViewModel();
ko.applyBindings(viewModel);
