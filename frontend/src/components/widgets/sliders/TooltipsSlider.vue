Base class for sliders.

Allows sliders with single or multivalues, + formatting of value and tooltips.

<template lang='pug'>
.container.slider-div
  .center-block(:id='sliderId', :style='{width: sliderWidth}')
</template>


<script>
import noUiSlider from 'nouislider'
import tools from '../../../tools'

export default {
  props: {
    sliderId: {
      type: String,
      default: function () { return tools.generateRandomID() }
    },
    min: {
      default: 0
    },
    max: {
      default: 10
    },
    valueFormat: {
      default: 0
    },
    tooltipsFormat: {
      default: 0
    },
    suffix: {
      default: ''
    },
    sliderWidth: {
      default: '120px'
    },
    step: {
      default: 1
    },
    tooltips: {
      default: true
    },
    connect: {
      default: true
    },
    value: {
      default: 1
    }
  },
  data: function () {
    var valueFormat = tools.nouiformatter(this.valueFormat, '')
    var format = tools.nouiformatter(this.tooltipsFormat, this.suffix)
    var ttformat = [format]
    if (this.value.length > 1) {
      for (var i = 1; i < this.value.length; i++) {
        ttformat.push(format)
      }
    }
    return {
      slider: null,
      innervalue: this.value,
      sliderConfig: {
        start: this.value,
        tooltips: ttformat,
        range: {
          'min': this.min,
          'max': this.max
        },
        step: parseFloat(this.step),
        format: valueFormat,
        connect: this.connect
      }
    }
  },
  mounted: function mounted () {
    var self = this
    this.slider = document.getElementById(this.sliderId)
    this.nouislider = noUiSlider.create(this.slider, this.sliderConfig)
    self.innervalue = this.nouislider.get()
    this.nouislider.on('change', function () {
      self.innervalue = tools.numerizeValue(this.get(), self.valueFormat)
    })
  },
  watch: {
    value: {
      handler: function (val) {
        this.nouislider.set(val)
        this.innervalue = val
      },
      deep: true
    },
    innervalue: {
      handler: function (val) {
        this.$emit('input', val)
      },
      deep: true
    }
  }
}
</script>

<style>

.noUi-horizontal {
	height: 6px;
  margin-top: 4px;
  margin-bottom: 18px;
}


.noUi-handle {
	width: 22px;
	height: 22px;
	left: -11px;
	top: -8px;
  border-radius: 15px
}

.noUi-handle:before,
.noUi-handle:after {
	content: "";
	display: none;
}

.noUi-connect {
	background: #0986d6;
	box-shadow: inset 0 0 3px rgba(51,51,51,0.45);
-webkit-transition: background 450ms;
	transition: background 450ms;
}

.slider-div {
  display: block;
  max-width: 90%;
  margin-top:10px;
  margin-bottom:10px;
}

.slid {
  float: right;
}

.noUi-tooltip {
  margin-bottom:-46px;
  border:none;
  background:none;
  font-size: 10px;
  font-family: 'Source Sans Pro';
  width: 50px
}

</style>
