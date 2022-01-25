y: pvFromw(wFromWetbulbDryBulb(wetbulbTemp, temp, self.totalPressure()), self.totalPressure()),

function pvFromw(w, totalPressure) {
    if (typeof w === "string") w = parseFloat(w);
    if (w < 0.000001) return 0;
    return totalPressure / (1 + 0.621945 / w);
}

function wFromWetbulbDryBulb(wetbulbTemp, temp, totalPressure) {
    var wsatWetBulb = satHumidRatioFromTempIp(wetbulbTemp, totalPressure);
    return ((1093 - 0.556 * wetbulbTemp) * wsatWetBulb - 0.24 * (temp - wetbulbTemp)) / (1093 + 0.444 * temp - wetbulbTemp);
}

function satHumidRatioFromTempIp(temp, totalPressure) {
    if (arguments.length !== 2) throw Error(`Not all parameters specified. temp: ${temp}; P: ${totalPressure}`);
    var satPress = satPressFromTempIp(temp);
    return (0.621945 * satPress) / (totalPressure - satPress);
}

// Saturation pressure in psi from temp in °F.
function satPressFromTempIp(temp) {
    const c8 = -1.0440397e4;
    const c9 = -1.129465e1;
    const c10 = -2.7022355e-2;
    const c11 = 1.289036e-5;
    const c12 = -2.4780681e-9;
    const c13 = 6.5459673;

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

