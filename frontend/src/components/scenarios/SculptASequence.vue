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
    el-checkbox(v-if='Object.keys(form.editedFeatures).length > 0' v-model='form.editFeatures') Edit features
    featureseditor(v-if='form.editedFeatures && form.editFeatures',
                      :sequence='form.sequence', v-model='form.editedFeatures')
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
import featureseditor from '../../components/widgets/FeaturesEditor/FeaturesEditor'

var bioparsers = require('bio-parsers')

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
        file: null,
        editFeatures: false,
        editedFeatures: {},
        sequence: ''
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
    digestionset,
    featureseditor
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
  watch: {
    'form.file.content': function (value) {
      var genbank = atob(value.split(',')[1])
      this.form.editedFeatures = {}
      var self = this
      bioparsers.genbankToJson(genbank, function (result) {
        if (result.length === 1) {
          var parsed = result[0].parsedSequence
          self.form.sequence = parsed.sequence
          parsed.features.map(function (feature, id) {
            console.log(feature)
            self.form.editedFeatures[id] = {
              id: id.toString(),
              start: feature.start,
              end: feature.end + 1,
              strand: feature.strand,
              label: feature.name,
              selected: false
            }
          })
        }
      })
    },
    'form.editedFeatures': {
      deep: true,
      handler: function (value) {
        console.log(value)
      }
    }
  },
  computed: {
    bars: function () {
      var data = this.queryStatus.polling.data
      if (!data) { return [] }
      return [
        {
          text: 'Constraint',
          index: data.constraint_index,
          total: data.constraint_total
        },
        {
          text: 'Objective',
          index: data.objective_index,
          total: data.objective_total
        },
        {
          text: 'Region',
          index: data.location_index,
          total: data.location_total
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