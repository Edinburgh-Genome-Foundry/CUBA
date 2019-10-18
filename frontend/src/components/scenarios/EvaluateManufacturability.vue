<template lang="pug">

.page
  h1  {{ infos.title }}
  img.icon.center-block(slot='title-img', :src='infos.icon')
  p.scenario-description.
    Submit a sequence(s), get plots of patterns impacting synthesis and assembly difficulty.
  web-links(:emailSubject="'[CUBA] Feedback on web app: ' + infos.title",
            tweetMessage="Find patterns that can impact your sequence's manufacturability:",
            :tweetUrl="'https://cuba.genomefoundry.org/' + infos.path")

  .form

    h4.formlabel Sequence(s) to analyze
    collapsible(title='Examples')
      file-example(filename='example_records.zip',
                    @input="function (e) { \
                      form.files = [e]; \
                    }",
                    fileHref='/static/file_examples/evaluate_manufacturability/example_records.zip',
                    imgSrc='/static/file_examples/generic_logos/linear_part_records.svg')
        p.
          Some example records to be checked for manufacturability-impacting patterns
    files-uploader(v-model='form.files', :multiple='true',
                      text="Drop multiple Genbank/Fasta (or click to select)")
    .files-number(v-if='form.files.length > 0')
      p {{ form.files.length }} file{{ form.files.length > 1 ? 's' : '' }}  selected


    backend-querier(:form='form', :backendUrl='infos.backendUrl',
                    :validateForm='validateForm', submitButtonText='Evaluate',
                    v-model='queryStatus')
    el-alert(v-if='queryStatus.requestError', :title="queryStatus.requestError",
       type="error", :closable="false")
  .results(v-if='!queryStatus.polling.inProgress && queryStatus.result.pdf_report')
    center
      download-button(text='Summary spreadsheet',
                     :filedata='queryStatus.result.spreadsheet')
      download-button(text='Annotated records',
                     :filedata='queryStatus.result.records')
      p(style='max-width: 500px;').
        Note: the feature labels indicate the constraints that are
        breached. For instance a "No BsaI" label means that there is indeed
        a BsaI site at this location, breaching the "No BsaI" constraint.
      iframe(:src="queryStatus.result.pdf_report.data" 
              style="width: 90%; max-width: 1000px; height:800px;" frameborder="0")
      

  powered-by(:softwareNames='infos.poweredby')
</template>

<script>
import learnmore from '../../components/widgets/LearnMore'

var infos = {
  title: 'Evaluate Manufacturability',
  navbarTitle: 'Evaluate Manufacturability',
  path: 'evaluate_manufacturability',
  description: '',
  backendUrl: 'start/evaluate_manufacturability',
  icon: require('../../assets/images/evaluate_manufacturability.svg'),
  poweredby: ['dnachisel']
}

export default {
  data () {
    return {
      form: {
        show_features: true,
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
    learnmore
  },
  infos: infos,
  methods: {
    handleSuccess (evt) {
      console.log(evt)
    },
    validateForm () {
      var errors = []
      if (this.form.files.length === 0) {
        errors.push('Provide at least one sequence.')
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
}

.el-select {
  width: 100%
}

.report-radio {
  text-align: center;
  margin: 0 auto;
  .el-radio {
    margin-bottom: 20px;
  }
}

.files-number p {
  margin-top: -10px;
  font-size: 0.8em;
}

.figure-preview {
  margin-bottom: 5em;
  img {
    width:100%;
  }
}
</style>
