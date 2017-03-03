<!-- src/components/Home.vue -->

<template lang="pug">
div
  h1 The EGF's Collection <br /> of Useful Biological Apps
  h2 Pick a scenario and start playing.
  scenariospanel
</template>

<script>
import scenariospanel from '../ScenariosPanel'
export default {
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
