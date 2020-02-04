<template lang="pug">
.page
  h1  {{ infos.title }}
  img.icon.center-block(slot='title-img', :src='infos.icon')
  p.scenario-description Extract features and save them as individual files.
  web-links(:emailSubject="'[CUBA] Feedback on web app: ' + infos.title",
            tweetMessage="Extract features and save them as individual files:",
            :tweetUrl="'https://cuba.genomefoundry.org/' + infos.path")

  .form

    h4.formlabel Upload sequence files
    files-uploader(v-model='form.files',
                   tip='Fasta or genbank, or zip', :multiple='true')
    p.inline Minimum sequence length of common parts (nt):
      el-input-number.inline(v-model="form.min_sequence_length", size="small",
                             :min=1, :max=10000)
    p.inline
      el-checkbox(v-model='form.direct_sense') Convert antisense features to direct sense in the exported files

    backend-querier(:form='form', :backendUrl='infos.backendUrl',
                    :validateForm='validateForm', submitButtonText='Extract features',
                    v-model='queryStatus')
    el-alert(v-show='queryStatus.requestError', :title="queryStatus.requestError",
       type="error", :closable="false")

  .results(v-if='!queryStatus.polling.inProgress && queryStatus.result')
    download-button(v-if='queryStatus.result.zip_file',
                    :filedata='queryStatus.result.zip_file',
                    text='Download Sequences')
    img.result_image(:src='queryStatus.result.figure_data')

  //- powered-by(:softwareNames='infos.poweredby')
</template>

<script>
import learnmore from '../../components/widgets/LearnMore'

var infos = {
  title: 'Feature Extractor',
  navbarTitle: 'Feature Extractor',
  path: 'feature_extractor',
  description: '',
  backendUrl: 'start/feature_extractor',
  icon: require('../../assets/images/extract_sequence_features.svg')
}

export default {
  data () {
    return {
      infos: infos,
      form: {
        files: [],
        min_sequence_length: 20,
        direct_sense: true
      },
      queryStatus: {
        polling: {},
        result: {},
        requestError: ''
      },
      goal_options: [
        {
          label: 'Extracted features',
          value: 'extracted_features'
        }
      ]
    }
  },
  components: {
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
