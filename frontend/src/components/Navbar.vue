<template lang="pug">
el-menu(mode='horizontal' @select="handleSelect")
  el-menu-item(index='home') Home
  el-submenu(index='2')
    template(slot='title') Scenarios
    el-menu-item(v-for='scenario in scenarios', :index="scenario.infos.path") {{scenario.infos.navbarTitle}}
  el-menu-item(index='about') About
</template>

<script>
import scenarios from '../pages/scenarios/scenarios.js'
console.log(scenarios.list[0].infos)
export default {
  data: () => ({ scenarios: scenarios.list }),
  methods: {
    handleSelect: function (key, keyPath) {
      this.$router.push(key)
    }
  }
}

</script>

<style scoped>
.el-menu {
  background-color: white;
  border-bottom: 2px solid #dddddd
}
</style>
