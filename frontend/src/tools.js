function round (value, precision) {
  var multiplier = Math.pow(10, precision || 0)
  return Math.round(value * multiplier) / multiplier
}

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
  round: round,
  nouiformatter: function (format, suffix) {
    return {
      to: function (value) {
        var val = (format.length ? format[parseInt(value)] : round(value, format))
        return val + suffix
      },
      from: function (value) { return value }
    }
  },
  extendArray: function (destination, source) {
    for (var property in source) {
      destination[property] = source[property]
    }
    return destination
  },
  valueTracker: {
    watch: {
      value: function (val) {
        console.log(val)
      }
    }
  },
  numerizeValue: numerizeValue
}
