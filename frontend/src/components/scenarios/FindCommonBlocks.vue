<template lang="pug">
.page
  h1  {{ infos.title }}
  img.icon.center-block(slot='title-img', :src='infos.icon')
  p.scenario-description Find sets of compatible overhangs for your assembly problem.
  web-links(:emailSubject="'[CUBA] Feedback on web app: ' + infos.title",
            tweetMessage="Find common regions between different DNA sequences:",
            :tweetUrl="'https://cuba.genomefoundry.org/' + infos.path")
  //- learnmore Bla bla bla

  .form

    h4.formlabel Upload sequence files
    collapsible(title='Examples')
      file-example(filename='Example sequences',
                   fileHref='/static/file_examples/find_common_blocks/sequences.zip',
                   @input='function (e) { form.files.push(e) }'
                   imgSrc='/static/file_examples/generic_logos/sequences_records.png')
    files-uploader(v-model='form.files',
                   tip='Fasta or genbank, or zip', :multiple='true')
    p.inline Minimal block size (in nucleotides):
      el-input-number.inline(v-model="form.min_block_size", size="small",
                             :min=10, :max=5000)
    p.inline
      span(style='margin-right: 15px') Block selection:
      el-radio(class='radio' v-model='form.block_selection' label='larger_first') Larger first
      el-radio(class='radio' v-model='form.block_selection' label='most_coverage_first') Most coverage first

    backend-querier(:form='form', :backendUrl='infos.backendUrl',
                    :validateForm='validateForm', submitButtonText='Find blocks',
                    v-model='queryStatus')
    el-alert(v-show='queryStatus.requestError', :title="queryStatus.requestError",
       type="error", :closable="false")

  .results(v-if='!queryStatus.polling.inProgress && queryStatus.result')
    download-button(v-if='queryStatus.result.zip_file',
                    :filedata='queryStatus.result.zip_file',
                    text='Download Sequences')
    img.result_image(:src='queryStatus.result.figure_data')


  powered-by(:softwareNames='infos.poweredby')
</template>

<script>
import learnmore from '../../components/widgets/LearnMore'

var infos = {
  title: 'Find Common DNA Blocks',
  navbarTitle: 'Find Common DNA Blocks',
  path: 'find-common-blocks',
  description: '',
  backendUrl: 'start/find_common_blocks',
  icon: require('../../assets/images/find_common_blocks.svg'),
  poweredby: ['geneblocks', 'dnafeaturesviewer']
}

export default {
  data () {
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
