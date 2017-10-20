<template lang="pug">
.page
  h1 {{infos.title}}
  img.icon.center-block(slot='title-img', :src='infos.icon' title='NOT a proto-nazi symbol !')
  p.center.
    Provide the results of an automated restriction digestion, and get a complete
    report on the results.


  .form
    h4.formlabel Analysis
    el-select(v-model='form.goal', placeholder='Select')
      el-option(v-for='item in goal_options', :label='item.label',
                :value='item.value', :key='item.value')



    h4.formlabel Constructs Sequences
      helper(help="The title of each Genbank file (or Fasta entry) must be the construct's name.")

    collapsible(title='Examples')
      file-example(filename='constructs_sequences.zip',
                   fileHref='/static/file_examples/analyze_digests/constructs_sequences.zip',
                   imgSrc='/static/file_examples/sequences_records.png')
        p.
          Collection of constructs assembled using random parts from the EMMA standard.
          Unzip the file and drag the genbank files into the file upload area.

    files-uploader(v-model='form.constructsSequences', help='Fasta or Genbank files')
    el-checkbox(v-model='form.circularSequences') Sequences are circular



    .animated.fadeIn(v-if="form.goal === 'validation'")

      h4.formlabel Constructs Map

        helper(help="A spreadsheet where the cells contents indicate the construct in each well")

      collapsible(title='Examples')
        file-example(filename='constructs_sequences.zip',
                     fileHref='/static/file_examples/analyze_digests/constructs_map.xlsx',
                     imgSrc='/static/file_examples/analyze_digests/constructs_map.png')
          p.
            Collection of constructs assembled using random parts from the EMMA standard.
            Unzip the file and drag the genbank files into the file upload area.

      files-uploader(v-model='form.constructsMap', help='Fasta or Genbank files', :multiple='false')



    el-checkbox(v-model='form.severalDigestionsPerClone') Some clones have several digestions in different wells
    .animated.fadeIn(v-if='form.severalDigestionsPerClone')
      h4.formlabel Clones Map
        helper(help="A spreadsheet where the cells contents indicate the ID of the clone in each well")
      files-uploader(v-model='form.clonesMap', help='Fasta or Genbank files', :multiple='false')



    h4.formlabel Digestions
    p
      el-checkbox(v-model='form.uniqueDigestion') The digestion is the same in all wells
    digestionselector.animated.fadeIn(v-if='form.uniqueDigestion' v-model='form.digestion')
    .animated.fadeIn(v-else)
      p Provide a spreadsheet map of all digestions.
        helper(help='If a digestion has several enzymes, separate them with a comma')
      files-uploader(v-model='form.digestionsMap', :multiple='false')

    h4.formlabel Fragment Analysis ZIP File
      helper(help='Provide a ZIP file of the AATI fragment analyzer (only machine supported at the moment)')
    collapsible(title='Examples')
      file-example(filename='fragment_analyzer_data.zip',
                   fileHref='/static/file_examples/analyze_digests/fragment_analyzer_data.zip',
                   imgSrc='/static/file_examples/analyze_digests/fragment_analyzer_data_screenshot.png')
        p.
          Fragment analysis data from the digestion (DraI+HindIII) of constructs
          assembled using random parts from the EMMA standard. The results are
          bad enough to be totally unpublishable so here they are.
    files-uploader(v-model='form.fragmentAnalysisArchive', :multiple='false')



    h4.formlabel Validation parameters
    p Tolerance
      helper(help='Distance at which two bands are considered different (in proportion of the ladder span)')
      el-slider(:min='0', :max='1', :step='0.05', :show-tooltip='true', v-model='form.tolerance')
    p Bands range
      helper(help='All bands outside this range will be ignored')
      el-slider(:min='0', :max='15000', :range='true', :step='10', :show-tooltip='true', v-model='form.bandsRange')
    el-checkbox(v-model='form.includeDigestionPlots') Include plots of cutting sites (slower)



    backend-querier(:form='form', :backendUrl='infos.backendUrl',
                    :validateForm='validateForm', submitButtonText='Generate report',
                    v-model='queryStatus')
    el-alert(v-if='queryStatus.requestError && !queryStatus.polling.inProgress',
             :title="queryStatus.requestError",
       type="error", :closable="false")
    .results(v-if='!queryStatus.polling.inProgress && queryStatus.result.zip_file')
      download-button(v-if='queryStatus.result.zip_file', text='Download report',
                      :filedata='queryStatus.result.zip_file')
  powered-by(:softwareNames='infos.poweredby')
</template>

<script>
import learnmore from '../../components/widgets/LearnMore'
import sequencesuploader from '../../components/widgets/SequencesUploader'
import digestionselector from '../../components/widgets/DigestionSelector'

var infos = {
  title: 'Analyze Digests',
  navbarTitle: 'Analyze Digests',
  path: 'analyze-digests',
  description: '',
  backendUrl: 'start/analyze_digests',
  icon: require('assets/images/analyze_digests.svg'),
  poweredby: ['bandwitch', 'bandwagon']
}

export default {
  data: function () {
    return {
      form: {
        constructsMap: null,
        constructsSequences: {},
        clonesMap: null,
        uniqueDigestion: true,
        severalDigestionsPerClone: false,
        digestion: [],
        digestionsMap: null,
        goal: 'validation',
        fragmentAnalysisArchive: null,
        circularSequences: true,
        tolerance: 0.1,
        bandsRange: [10, 15000],
        includeDigestionPlots: true
      },
      infos: infos,
      ladder_options: [
        {
          label: 'Ladder 100 bp - 4000 bp',
          value: '100_to_4k'
        },
        {
          label: 'Ladder 35 bp - 5 kbp (AATI)',
          value: '35_to_5k'
        },
        {
          label: 'Ladder 75 bp - 15 kbp (AATI)',
          value: '75_to_15k'
        }
      ],
      goal_options: [
        {
          label: 'Sequence validation',
          value: 'validation'
        },
        {
          label: 'Sequence identification',
          value: 'identification'
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
    sequencesuploader,
    learnmore,
    digestionselector
  },
  infos: infos,
  methods: {
    handleSuccess: function (evt) {
      console.log(evt)
    },
    validateForm: function () {
      var errors = []
      if (this.form.constructsSequences.length < 1) {
        errors.push('Provide constructs sequences')
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

.el-slider {
  display: inline-block;
  margin-left: 1em;
  margin-bottom: -0.7em;
  width: 50%;
}


.el-select {
  width: 100%
}

.el-input-number.inline {
  margin-bottom: -9px;
  margin-left: 10px;
  width: 130px;
}

p.loadData {
  font-size: 0.8em;
  margin-bottom: 0;
  cursor: pointer;
}

.bands-range {
  .el-slider {
    display: inline-block;
    width: 150px;
    margin-bottom: -12px;
    margin-left: 10px;
  }
}
</style>
