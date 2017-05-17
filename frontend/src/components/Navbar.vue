<template lang="pug">
el-menu(:mode="fullWidth > 500 ? 'horizontal' : 'vertical'" @select="handleSelect")
  .logo(v-if='fullWidth > 500')
    img(src='../assets/images/logo.png')
  el-menu-item(index='home') Home
  el-submenu(index='2')
    template(slot='title') Scenarios
    el-menu-item(v-for='scenario in scenarios', :index="scenario.infos.path") {{scenario.infos.navbarTitle}}
  el-menu-item(index='about') About
</template>

<script>
import scenarios from './scenarios/scenarios.js'
console.log(scenarios.list[0].infos)
export default {
  data: () => ({
    scenarios: scenarios.list,
    fullWidth: 0
  }),
  methods: {
    handleSelect: function (key, keyPath) {
      this.$router.push(key)
    },
    handleResize: function (event) {
      this.fullWidth = document.documentElement.clientWidth
    }
  },
  mounted: function () {
    window.addEventListener('resize', this.handleResize)
    this.handleResize()
  }
}

</script>

<style lang='scss' scoped>
.el-menu, .el-submenu {
  background-color: white;
  max-width: 90%;
}

.el-menu {
  border-bottom: 2px solid #dddddd;
}

.logo {
  // display: inline-block;
  float: left;
  margin-left: 20px;
  margin-right: 20px;
  height:60px;
}

.logo img {
  height:140%;
}
</style>
