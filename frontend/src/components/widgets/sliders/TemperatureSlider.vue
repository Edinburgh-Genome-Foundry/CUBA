A slider for defining a temperature (in Celsius) or a range of temperatures
(use a single value or a range [55, 65])

<script>
import slider from './TooltipsSlider'
export default {
  mixins: [slider],
  props: {
    min: {
      default: 50
    },
    max: {
      default: 80
    },
    step: {
      default: 0.5
    },
    value: {
      default: () => [55, 75]
    },
    valueFormat: {
      default: 1
    },
    tooltipsFormat: {
      default: 1
    },
    suffix: {
      default: ' °C'
    }
  }
}
</script>
