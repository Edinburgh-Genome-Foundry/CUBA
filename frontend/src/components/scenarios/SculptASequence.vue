<template lang="pug">
div
  h1  {{ infos.title }}
  img.icon.center-block(slot='title-img', :src='infos.icon')
  p.center Optimize a sequence with annotations representing constraints and objectives.
  learnmore Bla bla bla

  .form
    h4.formlabel Provide an annotated sequence
    filesuploader(v-model='form.file', text="Drop a single Genbank file (or click to select)",
                      :multiple='false', help='')
    span(v-if='form.file') Selected: <i>{{ form.file.name }}</i>
    backend-querier(v-if='!validateForm().length',
                    :form='form', :backendUrl='infos.backendUrl',
                    :validateForm='validateForm', submitButtonText='Sculpt',
                    v-model='queryStatus')
    progress-bars(:bars='bars')

  el-alert(v-if='queryStatus.requestError', :title="queryStatus.requestError",
     type="error", :closable="false")
  .results(v-if='!queryStatus.polling.inProgress')
    p.results-summary(v-if='queryStatus.result.summary',
                     v-html="queryStatus.result.summary")
    download-button(v-if='queryStatus.result.zip_file',
                    :filedata='queryStatus.result.zip_file')

  powered-by(:softwareNames='infos.poweredby')
</template>

<script>
import learnmore from '../../components/widgets/LearnMore'
import filesuploader from '../../components/widgets/FilesUploader'
import digestionset from '../../components/widgets/DigestionSelectorSet'

var infos = {
  title: 'Sculpt A Sequence',
  navbarTitle: 'Sculpt A Sequence',
  path: 'sculpt-a-sequence',
  description: '',
  backendUrl: 'start/sculpt_a_sequence',
  icon: require('assets/images/sculpt_a_sequence.svg'),
  poweredby: ['dnachisel', 'dnafeaturesviewer']
}

export default {
  data: function () {
    return {
      form: {
        file: null
      },
      infos: infos,
      queryStatus: {
        polling: {},
        result: {},
        requestError: ''
      }
    }
  },
  components: {
    filesuploader,
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
      if (!this.form.file) {
        errors.push('Provide at least one file.')
      }
      return errors
    }
  },
  computed: {
    bars: function () {
      var data = this.queryStatus.polling.data
      if (!data) { return [] }
      return [
        {
          text: 'Pass',
          index: data.iteration_ind,
          total: data.n_iterations
        },
        {
          text: 'Failing constraint',
          index: data.evaluation_ind,
          total: data.n_evaluations
        },
        {
          text: 'Breach',
          index: data.location_ind,
          total: data.n_locations
        }
      ]
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
