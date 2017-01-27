<!-- src/components/Home.vue -->

<template lang="pug">
div
  h1 Welcome to the App Demo !
  h2 Pick a scenario and start playing.
  scenariospanel
</template>

<script>
import scenariospanel from './scenarios/ScenariosPanel'
export default {
  components: {
    scenariospanel
  }
}
</script>

<style scoped>
h2 {margin-bottom: 60px}
</style>
