<template lang="pug">

div
  h1 {{infos.title}}
  img.icon.center-block(slot='title-img', :src='infos.icon')
  p.center Submit sequences and an enzymatic mixes, get migration predictions.

  .form
    h4.formlabel Ladder
    el-select(v-model='form.ladder', placeholder='Select')
      el-option(v-for='item in ladder_options', :label='item.label', :value='item.value')
    h4.formlabel Digestions
    digestionset(v-model='form.digestions')
    h4.formlabel Sequences
    sequencesuploader(v-model='form.files')
    el-checkbox(v-model='form.make_report') Generate report
    results-area(:form='form', :backendUrl='infos.backendUrl')
</template>

<script>
import learnmore from '../../components/widgets/LearnMore'
import sequencesuploader from '../../components/widgets/SequencesUploader'
import digestionset from '../../components/widgets/DigestionSelectorSet'

var infos = {
  title: 'Predict Digests',
  navbarTitle: 'Predict Digests',
  path: 'predict-digests',
  description: '',
  backendUrl: 'start/predict_digests',
  icon: require('assets/images/predict-icon.svg')
}

export default {
  data: function () {
    return {
      form: {
        ladder: '100-4k',
        digestions: [],
        make_report: false,
        files: []
      },
      infos: infos,
      ladder_options: [
        {
          label: 'Ladder 100 bp - 4000 bp',
          value: '100-4k'
        }
      ]
    }
  },
  components: {
    sequencesuploader,
    learnmore,
    digestionset
  },
  infos: infos,
  methods: {
    handleSuccess: function (evt) {
      console.log(evt)
    }
  }
}
</script>

<style scoped>

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

.el-checkbox {
  font-weight: normal;
}


.el-select {
  width: 100%
}
</style>
