A slider for defining a sequence length or a range of lengths
(use a single value or a range [1000, 3000])

<script>
import slider from './TooltipsSlider'
export default {
  mixins: [slider],
  props: {
    min: {
      default: 0
    },
    max: {
      default: 5000
    },
    step: {
      default: 10
    },
    value: {
      default: () => [1000, 3000]
    },
    valueFormat: {
      default: 0
    },
    tooltipsFormat: {
      default: 0
    },
    suffix: {
      default: ' bp'
    }
  }
}
</script>
