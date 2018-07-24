<template lang='pug'>
.progress-bars
  .bar(v-for="bar in order" v-if='bars[bar] && bars[bar].index <= bars[bar].total')
    div.bar-name {{bars[bar].title}} {{ bars[bar].index }} / {{ bars[bar].total }}
    el-progress(:percentage="100 * bars[bar].index / bars[bar].total", :show-text='false')
</template>


<script>
export default {
  props: {
    bars: {default: () => ({})},
    order: {default: () => ({})}
  }
}
</script>

<style lang='scss' scoped>
.progress-bars {
  .bar {
    margin-bottom: 1em;
    .bar-name {
      font-family: 'Inconsolata', Courier;
      margin-bottom: 0.3em;
    }
  }


  max-width: 500px;
  margin: 1em auto;
}
</style>
