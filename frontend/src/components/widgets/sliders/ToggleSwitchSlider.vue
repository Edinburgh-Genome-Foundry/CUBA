<script>
import slider from './TooltipsSlider'
export default {
  mixins: [slider],
  props: {
    min: {
      default: 0
    },
    max: {
      default: 1
    },
    step: {
      default: 1
    },
    value: {
      default: 1
    },
    valueFormat: {
      default: 0
    },
    tooltipsFormat: {
      default: function () { return ['disabled', 'enabled'] }
    },
    connect: {
      default: function () { return [true, false] }
    },
    sliderWidth: {
      default: '30px'
    }
  },
  watch: {
    value: function (val) {
      console.log('ts', val)
    }
  }
}
</script>
