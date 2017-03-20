<template lang="pug">

div
  h1 Digestion Pattern Predictor
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

    backend-querier(:form='form', :backendUrl='infos.backendUrl',
                    :validateForm='validateForm', submitButtonText='Predict patterns',
                    v-model='queryStatus')
    el-alert(v-if='queryStatus.requestError', :title="queryStatus.requestError",
       type="error", :closable="false")
    .results(v-if='!queryStatus.polling.inProgress')
      download-button(v-if='queryStatus.result.file',
                      :filedata='queryStatus.result.file')
      .results-summary(v-if='queryStatus.result.preview',
                       v-html="queryStatus.result.preview.html")
</template>

<script>
import learnmore from '../../components/widgets/LearnMore'
import sequencesuploader from '../../components/widgets/SequencesUploader'
import digestionset from '../../components/widgets/DigestionSelectorSet'

var infos = {
  title: 'Predict Digestions',
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
      ],
      queryStatus: {
        polling: {},
        result: {},
        requestError: ''
      }
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
    },
    validateForm: function () {
      var errors = []
      if (this.form.digestions.length === 0) {
        errors.push('Provide at least one digestion')
      }
      if (this.form.files.length === 0) {
        errors.push('Provide at least one sequence file')
      }
      return errors
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
