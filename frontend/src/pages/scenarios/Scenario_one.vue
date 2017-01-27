<template lang="pug">
scenario(:title='infos.title', v-model='fields', :backendUrl='infos.backendUrl')
  div(slot="head")
    img.title-img.center-block(slot='title-img', :src='infos.icon')
    | In this scenario we will learn how to make mammoths.
    el-collapse(slot='learnmore-paragr')
      el-collapse-item(title='Learn more', name='1')
        div
          | Bla bla bla.
          | Bla bla bla.

  div(slot='form')
    h4 NAME
    el-input(placeholder="Name you pet" v-model="fields.name")
    h4 Number of legs
    el-slider(v-model="fields.legs", :min='legsRange[0]', :max='legsRange[1]')
    h4 Number of eyes
    el-slider(v-model="fields.eyes", :min='eyesRange[0]', :max='eyesRange[1]')
</template>

<script>
import scenario from './Scenario'

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
      fields: {
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
    scenario
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
