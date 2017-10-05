<template lang="pug">
div
  h1  {{ infos.title }}
  img.icon.center-block(slot='title-img', :src='infos.icon')
  p.center Find sets of compatible overhangs for your assembly problem.
  web-links(:mailSubject="'[CUBA] Feedback on web app: ' + infos.title",
            tweetMessage="Find common regions between different DNA sequences:",
            :tweetUrl="'http://cuba.genomefoundry.org/' + infos.path")
  learnmore Bla bla bla

  .form

    h4.formlabel Upload sequence files
    filesuploader(v-model='form.files', text="Drop files (or click to select)",
                  help='Fasta or genbank. No file too large please :)', :multiple='true')
    p.inline Minimal block size in basepairs:
      el-input-number.inline(v-model="form.min_block_size", size="small",
                             :min=10, :max=5000)
    p.inline Block selection:
      el-radio(class='radio' v-model='form.block_selection' label='larger_first') Larger first
      el-radio(class='radio' v-model='form.block_selection' label='most_coverage_first') Most coverage first

    backend-querier(:form='form', :backendUrl='infos.backendUrl',
                    :validateForm='validateForm', submitButtonText='Evaluate',
                    v-model='queryStatus')
    el-alert(v-show='queryStatus.requestError', :title="queryStatus.requestError",
       type="error", :closable="false")

  .results(v-if='!queryStatus.polling.inProgress && queryStatus.result')
    download-button(v-if='queryStatus.result.zip_file',
                    :filedata='queryStatus.result.zip_file')
    img.result_image(:src='queryStatus.result.figure_data')


  powered-by(:softwareNames='infos.poweredby')
</template>

<script>
import learnmore from '../../components/widgets/LearnMore'
import filesuploader from '../../components/widgets/FilesUploader'

var infos = {
  title: 'Find common DNA blocks',
  navbarTitle: 'Find common DNA blocks',
  path: 'find-common-blocks',
  description: '',
  backendUrl: 'start/find_common_blocks',
  icon: require('assets/images/find_common_blocks.svg'),
  poweredby: ['geneblocks', 'dnafeaturesviewer']
}

export default {
  data: function () {
    return {
      infos: infos,
      form: {
        files: [],
        min_block_size: 60,
        block_selection: 'most_coverage_first'
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
    validateForm: function () {
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
