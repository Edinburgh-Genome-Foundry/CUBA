<template lang="pug">

.page
  h1  {{ infos.title }}
  img.icon.center-block(slot='title-img', :src='infos.icon')
  p.scenario-description.
    Submit a sequence(s), get plots of patterns impacting synthesis and assembly difficulty.
  web-links(:emailSubject="'[CUBA] Feedback on web app: ' + infos.title",
            tweetMessage="Find patterns that can impact your sequence's manufacturability:",
            :tweetUrl="'http://cuba.genomefoundry.org/' + infos.path")

  .form

    h4.formlabel Sequence(s) to analyze
    filesuploader(v-model='form.files', :multiple='true',
                      text="Drop multiple Genbank/Fasta (or click to select)")
    .files-number(v-if='form.files.length > 0')
      p {{ form.files.length }} file{{ form.files.length > 1 ? 's' : '' }}  selected
    .report-radio
      el-radio(v-model='form.report', class="radio", label='quick_view') Quick view
      el-radio(v-model='form.report', class="radio", label='pdf_report') PDF report


    backend-querier(:form='form', :backendUrl='infos.backendUrl',
                    :validateForm='validateForm', submitButtonText='Evaluate',
                    v-model='queryStatus')
    el-alert(v-if='queryStatus.requestError', :title="queryStatus.requestError",
       type="error", :closable="false")
    .results(v-if='!queryStatus.polling.inProgress')
      download-button(v-if='queryStatus.result.pdf_report',
                      :filedata='queryStatus.result.pdf_report')
      .results-summary(v-if='queryStatus.result.preview',
                       v-html="queryStatus.result.preview.html")
      .figures-preview(v-if='queryStatus.result.figures_data')
        .figure-preview(v-for='fig in queryStatus.result.figures_data')
          h4 {{fig.filename}}
          img(:src='fig.img_data')

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
      enzymes: ['BsaI', 'BsmBI', 'BbsI'],
      form: {
        report: 'quick_view',
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
