function round (value, precision) {
  var multiplier = Math.pow(10, precision || 0)
  return Math.round(value * multiplier) / multiplier
}

let f = (a, b) => [].concat(...a.map(a => b.map(b => [].concat(a, b))))
let cartesian = (a, b, ...c) => b ? cartesian(f(a, b), ...c) : a

function numerizeValue (value, rounding) {
  if (window.Array.isArray(value)) {
    var result = []
    value.forEach(function (val) {
      result.push(numerizeValue(val, rounding))
    })
    return result
  }
  var valuefloat = parseFloat(value)
  if (!isNaN(valuefloat)) {
    if (!isNaN(parseFloat(rounding))) {
      valuefloat = round(valuefloat, rounding)
    }
    return valuefloat
  } else {
    return value
  }
}

export default {
  generateRandomID: function () {
    return 'xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx'.replace(/[xy]/g, function (c) {
      var r = Math.random() * 16 | 0
      var v = c === 'x' ? r : (r & 0x3 | 0x8)
      return v.toString(16)
    })
  },
  numerizeValue: function (value, rounding) {
    if (window.Array.isArray(value)) {
      var result = []
      value.forEach(function (val) {
        result.push(numerizeValue(val, rounding))
      })
      return result
    }
    var valuefloat = parseFloat(value)
    if (!isNaN(valuefloat)) {
      if (!isNaN(parseFloat(rounding))) {
        valuefloat = round(valuefloat, rounding)
      }
      return valuefloat
    } else {
      return value
    }
  },
  cartesian: cartesian
}
