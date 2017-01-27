A slider for defining DNA per-basepair (by default, in dollars per basepair)

<script>
import slider from './TooltipsSlider'
export default {
  mixins: [slider],
  props: {
    min: {
      default: 0
    },
    max: {
      default: 1.0
    },
    step: {
      default: 0.01
    },
    value: {
      default: () => [0.20]
    },
    valueFormat: {
      default: 2
    },
    tooltipsFormat: {
      default: 2
    },
    suffix: {
      default: ' $/bp'
    }
  }
}
</script>
