<template lang="pug">
div
  h1 {{infos.title}}
  img.icon(slot='title-img', :src='infos.icon')
  p Lets make mammoths.
  learnmore
    p Bla bla bla

  .form
    h4 NAME
    el-input(placeholder="Name you pet" v-model="form.name")
    h4 Number of legs
    el-slider(v-model="form.legs", :min='legsRange[0]', :max='legsRange[1]')
    h4 Number of eyes
    el-slider(v-model="form.eyes", :min='eyesRange[0]', :max='eyesRange[1]')
    resultsarea(:form='form', :backendUrl='infos.backendUrl')
</template>

<script>
import resultsarea from './ResultsArea'
import learnmore from '../../components/widgets/LearnMore'

var infos = {
  title: 'Scenario 2',
  navbarTitle: 'Scenario 2',
  path: 'scenario-two',
  description: 'Design and order your mammoths',
  backendUrl: 'start/scenario_two',
  icon: require('assets/images/beaker-icon.svg')
}

export default {
  data: function () {
    return {
      form: {
        name: '',
        legs: 4,
        eyes: 2
      },
      infos: infos,
      eyesRange: [0, 10],
      legsRange: [0, 10]
    }
  },
  components: {
    resultsarea,
    learnmore
  },
  infos: infos
}
</script>

<style>
.el-collapse {
  border: none;
}

.el-collapse-item__header {
  padding-left: 0px;
}
</style>
