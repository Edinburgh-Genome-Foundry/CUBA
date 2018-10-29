<template lang="pug">
.page
  h1  {{ infos.title }}
  img.icon.center-block(slot='title-img', :src='infos.icon')
  p.scenario-description Find out what has been added/deleted/changed between two sequences.
  web-links(:emailSubject="'[CUBA] Feedback on web app: ' + infos.title",
            tweetMessage="Quick sequence difference visualization app:",
            :tweetUrl="'https://cuba.genomefoundry.org/' + infos.path")
  //- learnmore Bla bla bla

  .form
    h4.formlabel Sequence 1
    collapsible(title='Example')
      file-example(filename='sequence1.gb',
                   @input='function (e) { form.sequence1 = e }',
                   fileHref='/static/file_examples/compare_two_sequences/sequence1.gb')
    files-uploader(v-model='form.sequence1',
                   tip='Genbank/Fasta/Snapgene/txt/Zip',
                   :multiple='false')
    h4.formlabel Sequence 2
    collapsible(title='Example')
      file-example(filename='sequence2.gb',
                   @input='function (e) { form.sequence2 = e }'
                   fileHref='/static/file_examples/compare_two_sequences/sequence2.gb')
    files-uploader(v-model='form.sequence2',
                   tip="Genbank/Fasta/Snapgene/txt/Zip",
                   :multiple='false')

    el-form
      el-form-item(label='Figure width (in)')
        el-input-number(v-model='form.figure_width', :max='50', :min='4', :step='1', size='small')

    backend-querier(:form='form', :backendUrl='infos.backendUrl',
                    :validateForm='validateForm', submitButtonText='Evaluate',
                    v-model='queryStatus')
    el-alert(v-show='queryStatus.requestError', :title="queryStatus.requestError",
       type="error", :closable="false")

  .results(v-if='!queryStatus.polling.inProgress && queryStatus.result')
    center
      img.result_image(:src='queryStatus.result.figure_data')
    download-button(v-if='queryStatus.result.record',
                    text='Diff-annotated Genbank file',
                    :filedata='queryStatus.result.record')



  powered-by(:softwareNames='infos.poweredby')
</template>

<script>
import learnmore from '../../components/widgets/LearnMore'

var infos = {
  title: 'Compare Two Sequences',
  navbarTitle: 'Compare Two Sequences',
  path: 'compare-two-sequences',
  description: '',
  backendUrl: 'start/compare_two_sequences',
  icon: require('../../assets/images/compare_two_sequences.svg'),
  poweredby: ['geneblocks', 'dnafeaturesviewer']
}

export default {
  data () {
    return {
      infos: infos,
      form: {
        sequence1: null,
        sequence2: null,
        figure_width: 8
      },
      queryStatus: {
        polling: {},
        result: {},
        requestError: ''
      }
    }
  },
  components: {
    learnmore
  },
  infos: infos,
  methods: {
    validateForm () {
      var errors = []
      if (!this.form.sequence1.content) {
        errors.push('Provide a sequence 1 !')
      }
      if (!this.form.sequence2.content) {
        errors.push('Provide a sequence 2 !')
      }
      return errors
    }
  }
}
</script>

<style lang='scss' scoped>

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
  font-size: 16px !important;
}

.el-checkbox.inline {
  margin-left: 15px;
}


.el-slider.inline {
  display: inline-block;
  margin-left: 20px;
  margin-right: 20px;
  margin-bottom: -12px;
}

.el-input-number.inline {
  margin-bottom: -9px;
  margin-left: 10px;
  width: 130px;
}

.el-select {
  width: 100%
}

p.selected-overhangs {
  max-width: 90%;
  margin-left: 5%;

  span {
    font-family: 'Inconsolata', Courier;
    display: inline-block;
    margin: 0.2em 1em;
    padding: 3px;
    background-color: #eef3ff;
  }
}

.result_image {
  max-width: 80%;
  margin: 0 10%;
}
</style>
