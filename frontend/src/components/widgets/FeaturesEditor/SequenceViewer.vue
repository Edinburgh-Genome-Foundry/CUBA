<template lang='pug'>
.sequence-viewer
    svg(width='600', :viewBox="viewBox")
      text.nucl(v-for='nt, i in sequence.slice(window.start, window.end)',
               :x='i + window.start + 1 + 0.5', y='1', :key='i') {{nt}}
      text.ruler(v-for='nt, i in sequence.slice(window.start, window.end)',
                 :x='i +  window.start + 1 + 0.5', y='2', :key='i'
                 v-if='!((i +  window.start + 1) % 10)') {{i +  window.start + 1}}
</template>
<script>
import graphicfeature from './GraphicFeature'
export default {
  name: 'sequence-viewer',
  props: {
    sequence: {default: ''},
    window: {default: () => ({start: 0, end: 20})}
  },
  data () {
    return {
    }
  },
  components: {
    'graphic-feature': graphicfeature
  },
  computed: {
    viewBox () {
      return [
        this.window.start + 1,
        0,
        this.window.end - this.window.start,
        2.5
      ].join(' ')
    }
  }
}
</script>

<style lang='scss' scoped>

.sequence-viewer {

  width: 100%;
  svg {
    width: 100%;
    text {
      &.nucl {
        font-size: 1px;
      }
      &.ruler {
        font-size: 0.7px;
      }
      text-anchor: middle;
    }

    rect.focus {
      fill: rgba(0, 0, 0, 0.1);
    }
  }
}
</style>
