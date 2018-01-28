<template lang="pug">
.page
  el-alert(type='error' center :closable="false" show-icon) Currently under development
  h1  {{ infos.title }}
  img.icon.center-block(slot='title-img', :src='infos.icon')
  p.scenario-description Plot assembly features: highlight the important, discard the irrelevant.
  web-links(:emailSubject="'[CUBA] Feedback on web app: ' + infos.title",
            tweetMessage="Find common regions between different DNA sequences:",
            :tweetUrl="'http://cuba.genomefoundry.org/' + infos.path")
  //- learnmore Bla bla bla

  .form

    h4.formlabel Upload sequence files
    filesuploader(v-model='form.files', text="Drop files (or click to select)",
                  help='Fasta, Genbank, or Snapgene files. No file too large please :)', :multiple='true')
    h4.formlabel Plot style
    p.center
      el-radio(class='radio' v-model='form.linearity' label='linear') Linear
      el-radio(class='radio' v-model='form.linearity' label='circular') Circular
    p.inline.plot-size Width (inches)
      el-input-number.inline(v-model="form.plot_width", size="small",
                               :min='4', :max='16')
    p.inline.plot-size Height (inches)
      el-checkbox.inline(v-model='form.auto_height') auto
      el-input-number.inline(v-if='!form.auto_height', v-model="form.plot_height", size="small",
                                :min='4', :max='16')
    p
      el-checkbox(v-model='form.plot_full_sequence') Plot the full sequence's span
    p.inline(v-if='!form.plot_full_sequence') Plot index segment
      el-input-number.inline(v-model="form.plot_from_position", size="small",
                             :min='1', :max='form.plot_to_position-1')
      span &nbsp;to
      el-input-number.inline(v-model="form.plot_to_position", size="small",
                             :min='form.plot_from_position+1', :max='1000000')
    p
      el-checkbox(v-model='form.inline_labels') Allow inline feature labels.
    p
      el-checkbox(v-model='form.plot_sequence') Indicate sequence nucleotides
    p
      el-checkbox(v-model='form.plot_translation') Indicate amino acids

    h4.formlabel Ignore features where:
    el-button.center(icon='el-icon-plus' size='small') Add a rule

    h4.formlabel Custom feature styles
    p.center.inline Default feature color:
      span.options
        el-color-picker.inline(v-model="form.default_color", size="small")
    el-button.center(icon='el-icon-plus' size='small') Add a new style

    backend-querier(:form='form', :backendUrl='infos.backendUrl',
                    :validateForm='validateForm', submitButtonText='Plot',
                    v-model='queryStatus')
    el-alert(v-show='queryStatus.requestError', :title="queryStatus.requestError",
       type="error", :closable="false")

  .results(v-if='!queryStatus.polling.inProgress && queryStatus.result')
    download-button(v-if='queryStatus.result.pdf_file',
                    :filedata='queryStatus.result.pdf_file')


  powered-by(:softwareNames='infos.poweredby')
</template>

<script>
import learnmore from '../../components/widgets/LearnMore'
import filesuploader from '../../components/widgets/FilesUploader'

var infos = {
  title: 'Plot sequence features',
  navbarTitle: 'Plot sequence features',
  path: 'plot-common-features',
  description: '',
  backendUrl: 'start/plot_common_features',
  icon: require('assets/images/plot_sequence_features.svg'),
  poweredby: ['dnafeaturesviewer']
}

export default {
  data () {
    return {
      infos: infos,
      form: {
        files: [],
        linearity: 'linear',
        default_color: '#BABAEB',
        plot_width: 12,
        plot_height: 6,
        auto_height: true,
        inline_labels: false,
        plot_full_sequence: true,
        plot_from_position: 1,
        plot_to_position: 10000
      },
      queryStatus: {
        polling: {},
        result: {},
        requestError: ''
      },
      goal_options: [
        {
          label: 'A collection of compatible overhangs',
          value: 'overhangs_set'
        },
        {
          label: 'A sequence decomposition, with compatible overhangs',
          value: 'sequence_decomposition'
        }
      ]
    }
  },
  components: {
    filesuploader,
    learnmore
  },
  infos: infos,
  methods: {
    validateForm () {
      var errors = []
      if (!this.form.files.length) {
        errors.push('Provide files !')
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

.form /deep/ {
  .el-color-picker__trigger {
    margin-bottom: -9px;
  }
  .el-checkbox__label {
    font-weight: normal;
    color: #2c3e50;
    font-size: 16px;
  }
}

.title-img {
  height:80px;
  margin-top: -20px;
  margin-bottom: 20px;
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

.plot-size {
  .el-input-number {
    width: 100px;
  }
}

.el-select {
  width: 100%
}

.result_image {
  max-width: 80%;
  margin: 0 10%;
}
</style>
