<!-- src/components/Home.vue -->

<template lang="pug">
div
  h1 The EGF's Collection <br /> of Useful Biological Apps
  h2 Pick a scenario below.
  .scenario-category(v-for='category in scenarios', :key='category')
    h3 {{category.category}}
    scenariospanel(:scenarios='category.scenarios')
</template>

<script>
import scenariospanel from '../widgets/ScenariosPanel'
import scenarios from '../scenarios/scenarios.js'
export default {
  data: function () {
    return {
      scenarios: scenarios.list
    }
  },
  components: {
    scenariospanel
  }
}
</script>

<style scoped>
h2 {margin-bottom: 60px}
p {
  font-size: 18px;
}
</style>
