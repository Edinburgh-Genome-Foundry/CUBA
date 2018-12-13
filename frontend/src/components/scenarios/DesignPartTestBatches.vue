<template lang="pug">
.page
  h1  {{ infos.title }}
  img.icon.center-block(slot='title-img', :src='infos.icon')
  p.scenario-description.
    If you are going to do many (combinatorial) assemblies in parralel, this app
    can help you select a few assemblies to run first to test no part is
    corrupted.
  web-links(:emailSubject="'[CUBA] Feedback on web app: ' + infos.title",
            tweetMessage="Sequenticons are human-friendly, visual DNA sequence identifiers.",
            :tweetUrl="'https://cuba.genomefoundry.org/' + infos.path")
  //- learnmore Bla bla bla

  .form
    
    el-select(v-model='form.input_format')
      el-option(value='list' label='I will provide a list of assemblies')
      el-option(value='combinatorial' label='I will provide a combinatorial design')
    p Max. expected number of saboteur parts: &nbsp; &nbsp;
      el-input-number(v-model='form.max_saboteurs', :min='1', :max='10' size='mini')
    h4.formlabel {{ form.input_format === 'list' ? 'Assembly list' : 'Combinatorial design'}}
    collapsible(title='Examples')
      file-example(:filename='`${form.input_format}_assembly_plan.csv`',
                   :fileHref='`/static/file_examples/design_part_test_batches/${form.input_format}_assembly_plan.csv`',
                   @input='function (e) { form.input_file = e }'
                   imgSrc='/static/file_examples/generic_logos/spreadsheet.svg')
    files-uploader(v-model='form.input_file', text="Drop files (or click to select)",
                  help='CSV files only :)', :multiple='false')

    backend-querier(:form='form', :backendUrl='infos.backendUrl',
                    :validateForm='validateForm', submitButtonText='Design a test batch',
                    v-model='queryStatus')
    el-alert(v-show='queryStatus.requestError', :title="queryStatus.requestError",
       type="error", :closable="false")

  .results(v-if='!queryStatus.polling.inProgress && Object.keys(queryStatus.result).length')
    p(align='center').
      A batch with #[b {{ queryStatus.result.n_selected }}] selected assemblies was designed.
    download-button(v-if='queryStatus.result.report',
                    text='Download the report',
                    :filedata='queryStatus.result.report')


  powered-by(:softwareNames='infos.poweredby')
</template>

<script>
import learnmore from '../../components/widgets/LearnMore'

var infos = {
  title: 'Design Part Test Batches',
  navbarTitle: 'Design Part Test Batches',
  path: 'design_part_test_batches',
  description: '',
  backendUrl: 'start/design_part_test_batches',
  icon: require('../../assets/images/design_part_test_batches.svg'),
  poweredby: ['saboteurs']
}

export default {
  data () {
    return {
      infos: infos,
      form: {
        input_format: 'list',
        input_file: null,
        max_saboteurs: 1
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
      if (!this.form.input_file) {
        errors.push('Provide files !')
      }
      return errors
    }
  }
}
</script>

<style lang='scss' scoped>

.form {
  margin: 50px auto;
  max-width: 500px;
}

h4.formlabel {
  text-align: center;
  text-transform: uppercase;
  margin-top: 40px
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
