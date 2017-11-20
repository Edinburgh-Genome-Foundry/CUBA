<template lang='pug'>
.features-viewer
    svg(width='600', :height='30*(maxLevel + 2)', :viewBox="viewBox",
        preserveAspectRatio="none",
        @mousedown="startDrag", @touchstart="startDrag",
        @mousemove="onDrag", @touchmove="onDrag",
        @mouseup="stopDrag", @touchend="stopDrag",
        @mouseleave="stopDrag")
      path(:d="'M 0 0 L ' + sequenceLength + ' 0'" stroke='black', stroke-width='0.05px')
      graphic-feature(v-for='featureData, i in featuresData', :featureData='featureData',
                      :isSelected='featureData.selected',
                      :key='i', @click="$emit('select', featureData.id)",
                      @mouseover="function (v) {$emit('mouseover', v)}")
</template>
<script>
import graphicfeature from './GraphicFeature'
export default {
  name: 'features-viewer',
  props: {
    featuresData: {default: () => ([])},
    sequenceLength: {default: 20},
    window: {default: () => ({start: 0, end: 20})}
  },
  data () {
    return {
      dragging: false,
      originalX: null
    }
  },
  components: {
    'graphic-feature': graphicfeature
  },
  computed: {
    maxLevel () {
      if (this.featuresData.length === 0) {
        return 1
      }
      var result = Math.max.apply(1, this.featuresData.map((f) => f.level))
      return result
    },
    viewBox () {
      return [
        this.window.start,
        -0.5,
        this.window.end - this.window.start,
        this.maxLevel + 2
      ].join(' ')
    }
  },
  methods: {
    startDrag (evt) {
      // var e = evt.target
      var dim = this.$el.getBoundingClientRect()
      this.originalX = Math.round(Number((this.window.end - this.window.start) * (evt.clientX - dim.left) / (dim.right - dim.left)))
      this.dragging = true
    },
    onDrag (evt) {
      if (this.dragging) {
        var dim = this.$el.getBoundingClientRect()
        // var windowCenter = Math.round(Number(this.window * (evt.clientX - dim.left) / (dim.right - dim.left)))
        this.newX = Math.round(Number((this.window.end - this.window.start) * (evt.clientX - dim.left) / (dim.right - dim.left)))
        this.$emit('drag', this.originalX - this.newX)
        // console.log()
        this.originalX = this.newX
      }
    },
    stopDrag (evt) {
      this.dragging = false
    }
  }
}
</script>
<!-- Add "scoped" attribute to limit CSS to this component only -->
<style lang='scss' scoped>
.features-viewer {
  width: 100%;
  svg {
    width: 100%;
  }
}
</style>
