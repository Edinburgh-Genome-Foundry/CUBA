<template lang="pug">
el-menu(:mode="fullWidth > 610 ? 'horizontal' : 'vertical'" @select="handleSelect",
        :unique-opened='true')
  .logo(v-if='fullWidth > 610')
    router-link(to='home')
      img(src='../assets/images/cuba-title.png')
  el-menu-item(index='home') Home
  el-menu-item(index='about') About
  el-submenu(index='2')
    template(slot='title') Scenarios
    //- el-submenu(:index="'2-' + (index + 1)" v-for='(category, index) in scenarios')
    //-     span(slot='title') {{category.category}}
    div(v-for='category in scenarios')
      el-menu-item(v-for='(scenario, subindex) in category.scenarios',
                   :index="scenario.infos.path",
                   :key='scenario.infos.path') {{scenario.infos.navbarTitle}}
</template>

<script>
import scenarios from './scenarios/scenarios.js'
export default {
  data () {
    return {
      scenarios: scenarios.list,
      fullWidth: 0
    }
  },
  methods: {
    handleSelect (key, keyPath) {
      this.$router.push(key)
    },
    handleResize (event) {
      this.fullWidth = document.documentElement.clientWidth
    }
  },
  mounted () {
    window.addEventListener('resize', this.handleResize)
    this.handleResize()
  }
}
</script>

<style lang='scss'>
.el-submenu {
  height: 70px;
}

.logo {
  // display: inline-block;
  float: left;
  margin-left: 20px;
  margin-right: 20px;
  height:70px;
}


.logo img {
  height:140%;
}
</style>
