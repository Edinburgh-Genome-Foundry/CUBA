<template lang="pug">

.page
  h1 Digestion Pattern Predictor
  web-links(:emailSubject="'[CUBA] Feedback on web app: ' + infos.title",
            tweetMessage="Online restriction digest simulator:",
            :tweetUrl="'https://cuba.genomefoundry.org/' + infos.path")
  img.icon.center-block(slot='title-img', :src='infos.icon')
  p.scenario-description Submit sequences and an enzymatic mixes, get migration predictions.


  .form
    h4.formlabel Ladder
    ladderselector(v-model='form.ladder')

    h4.formlabel Digestions
    .example(@click="form.digestions = [['PvuI'], ['EcoRI', 'XbaI'], ['ApaLI', 'HindIII']]",
             style='cursor: pointer; color: grey; margin: 0.5em;') Click me to load an example
    digestionset(v-model='form.digestions')
    h4.formlabel Sequences
    collapsible(title='Examples')
      file-example(filename='emma_constructs.zip',
                  @input='function (e) {e.circularity = true; form.files = [e]}',
                  fileHref='/static/file_examples/predict_digestions/emma_constructs.zip',
                  imgSrc='/static/file_examples/generic_logos/part.svg')
    files-uploader(v-model='form.files', help='Fasta or Genbank files')
    el-checkbox(v-model='form.circular_sequences') Sequences are circular
    p
      el-checkbox(v-model='form.use_file_names_as_ids') Use file names as sequence IDs
    p
      el-checkbox(v-model='form.make_cuts_position_report') Generate report
    p
      el-checkbox(v-model='form.show_band_sizes') Show band sizes in plot.
    p
      el-checkbox(v-model='form.use_ordering_list') Provide a sequences ordering list
    div(v-if='form.use_ordering_list')
      h4.formlabel Sequences Plotting Order
      el-input(type='textarea', v-model='form.ordering_list', :rows=3,
               placeholder='Write sequence names in new lines or separated by commas')

    backend-querier(:form='form', :backendUrl='infos.backendUrl',
                    :validateForm='validateForm', submitButtonText='Predict patterns',
                    v-model='queryStatus')
    el-alert(v-if='!queryStatus.polling.inProgress && queryStatus.requestError', :title="queryStatus.requestError",
       type="error", :closable="false")
  .results(v-if='!queryStatus.polling.inProgress && queryStatus.result.figure_data')
    center
      img(:src='queryStatus.result.figure_data')
      iframe(:src="queryStatus.result.pdf_data" 
              style="width: 90%; max-width: 1000px; height:800px;"
              frameborder="0")
  powered-by(:softwareNames='infos.poweredby')
</template>

<script>
import learnmore from '../../components/widgets/LearnMore'
import sequencesuploader from '../../components/widgets/SequencesUploader'
import digestionset from '../../components/widgets/DigestionSelectorSet'
import ladderselector from '../../components/widgets/LadderSelector'

var infos = {
  title: 'Predict Digestions',
  navbarTitle: 'Predict Digestions',
  path: 'predict-digestions',
  description: '',
  backendUrl: 'start/predict_digestions',
  icon: require('../../assets/images/predict-icon.svg'),
  poweredby: ['bandwagon']
}

export default {
  data () {
    return {
      form: {
        ladder: '100_to_4k',
        circular_sequences: true,
        digestions: [],
        make_cuts_position_report: false,
        files: [],
        show_band_sizes: false,
        use_ordering_list: false,
        use_file_names_as_ids: true,
        ordering_list: ''
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
    sequencesuploader,
    learnmore,
    digestionset,
    ladderselector
  },
  infos: infos,
  methods: {
    handleSuccess (evt) {
      console.log(evt)
    },
    validateForm () {
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

<style lang='css' scoped>

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
