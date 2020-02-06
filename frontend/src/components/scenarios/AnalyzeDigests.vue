<template lang="pug">
.page
  h1 {{infos.title}}
  web-links(:emailSubject="'[CUBA] Feedback on web app: ' + infos.title",
            tweetMessage="Convert sequences formats to Genbank or Fasta",
            :tweetUrl="'https://cuba.genomefoundry.org/' + infos.path")
  img.icon.center-block(slot='title-img', :src='infos.icon')
  p.scenario-description.
    Provide the results of an automated restriction digestion, and get a complete
    report comparing your results to the expected patterns.


  .form
    h4.formlabel Analysis
    el-select(v-model='form.goal', placeholder='Select')
      el-option(v-for='item in goal_options', :label='item.label',
                :value='item.value', :key='item.value')

    div(v-if="form.goal === 'validation'", style='margin-top:1em;')
      el-radio(v-model='form.subanalysis', class="radio", label='standard') Standard analysis
      el-radio(v-model='form.subanalysis', class="radio", label='partial_digests')  Partial digests analysis



    h4.formlabel Constructs Sequences
      helper.
        The title of each Genbank file (or Fasta entry) must be the construct's name.

    collapsible(title='Examples')
      file-example(filename='constructs_sequences.zip',
                   fileHref='/static/file_examples/analyze_digests/constructs_sequences.zip',
                   @input='function (e) { form.constructsSequences.push(e) }',
                   imgSrc='/static/file_examples/sequences_records.png')
        p.
          Collection of constructs assembled using random parts from the EMMA standard.
          Unzip the file and drag the genbank files into the file upload area.

    files-uploader(v-model='form.constructsSequences', help='Fasta/Genbank/Zip files')
    p
      el-select(v-model='form.topology' size='small')
        el-option(value='circular' label='All sequences are circular')
        el-option(value='linear' label='All sequences are linear')
        el-option(value='default_to_circular' label='Autodetect each sequence\'s topology  (default to circular)')
        el-option(value='default_to_linear' label='Autodetect each sequence\'s topology (default to linear)')
    



    .animated.fadeIn(v-if="form.goal === 'validation'")

      h4.formlabel Constructs Map

        helper.
          A spreadsheet where the cells contents indicate the construct in each well

      collapsible(title='Examples')
        file-example(filename='constructs_map.xlsx',
                     @input='function (e) { form.constructsMap = e }',
                     fileHref='/static/file_examples/analyze_digests/constructs_map.xlsx',
                     imgSrc='/static/file_examples/analyze_digests/constructs_map.png')
          p.
            Collection of constructs assembled using random parts from the EMMA standard.
            Unzip the file and drag the genbank files into the file upload area.

      files-uploader(v-model='form.constructsMap', help='CSV or Excel file', :multiple='false')



    el-checkbox(v-model='form.severalDigestionsPerClone') Some clones have several digestions in different wells
    .animated.fadeIn(v-if='form.severalDigestionsPerClone')
      h4.formlabel Clones Map
        helper.
          A spreadsheet where the cells contents indicate the ID of the clone in each well
      files-uploader(v-model='form.clonesMap', help='CSV or Excel file', :multiple='false')



    h4.formlabel Digestions
    p
      el-checkbox(v-model='form.uniqueDigestion') The digestion is the same in all wells
    digestionselector.animated.fadeIn(v-if='form.uniqueDigestion' v-model='form.digestion')
    .animated.fadeIn(v-else)
      p Provide a spreadsheet map of all digestions.
        helper.
          If a digestion has several enzymes, separate them with a comma
      files-uploader(v-model='form.digestionsMap', :multiple='false')

    h4.formlabel Fragment Analysis ZIP File
      helper.
        Provide a ZIP file of the AATI fragment analyzer (only machine supported at the moment)
    collapsible(title='Examples')
      file-example(filename='fragment_analyzer_data.zip',
                   @input='function (e) { form.fragmentAnalysisArchive = e }',
                   fileHref='/static/file_examples/analyze_digests/fragment_analyzer_data.zip',
                   imgSrc='/static/file_examples/analyze_digests/fragment_analyzer_data_screenshot.png')
        p.
          Fragment analysis data from the digestion (DraI+HindIII) of constructs
          assembled using random parts from the EMMA standard. The results are
          bad enough to be totally unpublishable so here they are.
    files-uploader(v-model='form.fragmentAnalysisArchive', :multiple='false')



    h4.formlabel Validation parameters
    p Tolerance
      helper.
        Distance at which two bands are considered different (in proportion of the ladder span)
      el-input-number.inline(:min='0', :max='0.5', :step='0.05', :show-tooltip='true', v-model='form.tolerance' size='small')
    p Bands range
      helper.
        All bands outside this range will be ignored
      span From
      el-input-number.inline(:min='0', :max='form.bandsRange[1]',
                             v-model='form.bandsRange[0]', size='small', :step='10')
      span to
      el-input-number.inline(:min='form.bandsRange[0]', :max='15000',
                             v-model='form.bandsRange[1]', size='small', :step='10')
    p Ignore bands under
      el-input-number.inline(:min='0', :max='5000',
                           v-model='form.ignoreBandsUnder', size='small', :step='10')
      span bp
    el-checkbox(v-model='form.includeDigestionPlots') Include plots of cutting sites (slower)



    backend-querier(:form='form', :backendUrl='infos.backendUrl',
                    :validateForm='validateForm', submitButtonText='Generate report',
                    v-model='queryStatus')
    el-alert(v-if='queryStatus.requestError && !queryStatus.polling.inProgress',
             :title="queryStatus.requestError",
       type="error", :closable="false")
    .results(v-if='!queryStatus.polling.inProgress && (queryStatus.result.zip_file || queryStatus.result.message || queryStatus.result.pdf_file)')
      download-button(v-if='queryStatus.result.zip_file', text='Download report',
                      :filedata='queryStatus.result.zip_file')
      p(v-if='queryStatus.result.message', v-html='queryStatus.result.message')
      download-button(v-if='queryStatus.result.pdf_file', text='Download report',
                      :filedata='queryStatus.result.pdf_file')
      .center(v-if='queryStatus.result.figure_data')
        img(:src='queryStatus.result.figure_data')
    .unknown-constructs(v-if='queryStatus.result.unknown_constructs')
      el-alert(title="Couldn't find a corresponding record for some constructs in the plate."
               type="error")
      ul
        li(v-for='data, cst in queryStatus.result.unknown_constructs').
          #[b {{ cst }}] (well {{ data.well }}). <br/> Did you mean: {{ data.suggestions.join(', ')}}
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
  icon: require('../../assets/images/analyze_digests.svg'),
  poweredby: ['bandwitch', 'bandwagon']
}

export default {
  data () {
    return {
      form: {
        constructsMap: null,
        constructsSequences: [],
        clonesMap: null,
        uniqueDigestion: true,
        severalDigestionsPerClone: false,
        digestion: [],
        digestionsMap: null,
        goal: 'validation',
        fragmentAnalysisArchive: null,
        tolerance: 0.05,
        bandsRange: [10, 15000],
        includeDigestionPlots: true,
        subanalysis: 'standard',
        ignoreBandsUnder: 0,
        topology: 'default_to_circular'
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
    handleSuccess (evt) {
      console.log(evt)
    },
    validateForm () {
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
  margin-right: 10px;
  width: 120px;
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
