<template lang='pug'>
g.graphic-feature
  path(v-if="shape === 'arrow'", :d='path', @click="$emit('click')" v-bind:class="{ active: isSelected }",
       :fill='featureData.js_color',
        @mouseover="$emit('mouseover', featureData.name)",
        @mouseleave="$emit('mouseover', null)")
  rect(v-if="shape === 'rect'", :x='featureData.start', :y='featureData.level - 0.5 + 0.3',
       :width='featureData.end - featureData.start', :height='0.6', :rx='20', :ry='20',
       @click="$emit('click')" v-bind:class="{ active: isSelected }", :fill='featureData.js_color',
       @mouseover="$emit('mouseover', featureData.id)", @mouseleave="$emit('mouseover', null)")
</template>
<script>
export default {
  name: 'graphic-feature',
  props: {
    featureData: {default: () => ({})},
    isSelected: {default: false},
    showLabel: {default: false},
    shape: {default: 'arrow'}
  },
  computed: {
    path () {
      console.log(this.useShape)
      var points = []

      var middle = this.featureData.level

      var start = this.featureData.start - 1
      var end = this.featureData.end
      var top = this.featureData.level + 0.3
      var bottom = this.featureData.level - 0.3

      if (this.featureData.strand === 0) {
        points = [
          [start, top],
          [end, top],
          [end, bottom],
          [start, bottom]
        ]
      } else if (this.featureData.strand === 1) {
        points = [
          [start, top],
          [end - 1, top],
          [end, middle],
          [end - 1, bottom],
          [start, bottom]
        ]
      } else if (this.featureData.strand === -1) {
        points = [
          [end, top],
          [start + 1, top],
          [start, middle],
          [start + 1, bottom],
          [end, bottom]
        ]
      }
      var result = points.map(c => ['L', c[0], c[1]].join(' ')).join(' ') + ' Z'
      result = 'M' + result.slice(1)
      return result
    }
  }
}
</script>
<!-- Add 'scoped' attribute to limit CSS to this component only -->
<style lang='scss' scoped>
.graphic-feature {
  path, rect {
    &.active {
      stroke-width: 0.12px;
    }
    stroke: black;
    stroke-width: 0.05px;
    // fill: #aaaaff;
  }
}
</style>
