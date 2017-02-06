<template lang="pug">
div
  h1 {{infos.title}}
  img.title-img.center-block(slot='title-img', :src='infos.icon')
  | Lets make dinosaurs.
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
  title: 'Scenario 1',
  navbarTitle: 'Scenario 1',
  path: 'scenario-one',
  description: 'Design and order your dinosaurs',
  backendUrl: 'start/scenario_one',
  icon: require('assets/images/pencil-icon.svg')
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

h4.formlabel {
  text-align: center;
  text-transform: uppercase;
  margin-top: 40px
}

.form {
  margin: 50px auto;
  max-width: 500px;
}

.title-img {
  height:80px;
  margin-top: -20px;
  margin-bottom: 20px;

}
</style>
