<template lang='pug'>
.features-minimap
    svg(width='600', :height='15 * (maxLevel + 2)',
        preserveAspectRatio="none",
        :viewBox="viewBox",
        @mousedown="startDrag", @touchstart="startDrag",
        @mousemove="onDrag", @touchmove="onDrag",
        @mouseup="stopDrag", @touchend="stopDrag",
        @mouseleave="stopDrag")
      path(:d="'M 0 0 L ' + sequenceLength + ' 0'" stroke='black', stroke-width='0.05px')
      rect.focus(:x='window.start', :y='-1', :width='window.end - window.start', :height='maxLevel + 2')
      //- el-tooltip(v-for='featureData, i in featuresData', :content='featureData.label', :key='i')
      graphic-feature(v-for='featureData, i in featuresData', :key='i',
                      :featureData='featureData', :isSelected='featureData.selected',
                      @click="$emit('select', featureData.id)", :shape="'rect'",
                      @mouseover="function (v) {$emit('mouseover', v)}")

</template>
<script>
import graphicfeature from './GraphicFeature'
export default {
  name: 'minimap',
  props: {
    featuresData: {default: () => ([])},
    sequenceLength: {default: 20},
    window: {default: () => ({start: 0, end: 20})}
  },
  data () {
    return {
      dragging: false
    }
  },
  components: {
    'graphic-feature': graphicfeature
  },
  computed: {
    maxLevel () {
      console.log(this.featuresData)
      if (this.featuresData.length === 0) {
        return 1
      }
      var result = Math.max.apply(1, this.featuresData.map((f) => f.level))
      return result
    },
    viewBox () {
      return [0, -1, this.sequenceLength, this.maxLevel + 2].join(' ')
    }
  },
  methods: {
    startDrag (evt) {
      // var e = evt.target
      var dim = this.$el.getBoundingClientRect()
      var windowCenter = Math.round(Number(this.sequenceLength * (evt.clientX - dim.left) / (dim.right - dim.left)))
      this.dragging = true
      this.$emit('drag', windowCenter)
    },
    onDrag (evt) {
      if (this.dragging) {
        var dim = this.$el.getBoundingClientRect()
        var windowCenter = Math.round(Number(this.sequenceLength * (evt.clientX - dim.left) / (dim.right - dim.left)))
        this.$emit('drag', windowCenter)
      }
    },
    stopDrag (evt) {
      this.dragging = false
    },
    featureMouseover (evt) {
      console.log(evt)
    }
  }
}
</script>
<!-- Add "scoped" attribute to limit CSS to this component only -->
<style lang='scss' scoped>
.features-minimap {
  width: 100%;
  svg {
    width: 100%;

    rect.focus {
      fill: rgba(0, 0, 0, 0.1);
    }
  }
}
</style>
