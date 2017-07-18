<template lang='pug'>
.progress-bars
  p(v-for="bar in bars" v-if='bar.index < bar.total')
    span {{bar.text}} {{ bar.index }} / {{ bar.total }}
    el-progress(:percentage="100 * bar.index / bar.total", :show-text='false')
</template>


<script>
export default {
  props: {
    bars: {default: () => ([])}
  }
}
</script>

<style lang='scss' scoped>
.progress-bars {
  font-family: 'Inconsolata', Courier;
}
</style>
