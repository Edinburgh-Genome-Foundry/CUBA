<script>
import slider from './TooltipsSlider'
export default {
  mixins: [slider],
  props: {
    min: {
      default: 0
    },
    max: {
      default: 60
    },
    step: {
      default: 1
    },
    value: {
      default: 30
    },
    valueFormat: {
      default: 0
    },
    tooltipsFormat: {
      default: 0
    },
    suffix: {
      default: ' days'
    }
  }
}
</script>
