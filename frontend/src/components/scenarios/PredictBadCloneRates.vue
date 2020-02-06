<template lang="pug">

.page
  h1 Predict Bad Clone Rates
  web-links(:emailSubject="'[CUBA] Feedback on web app: ' + infos.title",
            tweetMessage="Online restriction digest simulator:",
            :tweetUrl="'https://cuba.genomefoundry.org/' + infos.path")
  img.icon.center-block(slot='title-img', :src='infos.icon')
  p.scenario-description.
    Find out how much unwanted parts cross-talk there will be in your next assemblies,
    and which proportion of clones will have the right constructs. Warning: this is
    yet very experimental.
  el-alert(title="This is a very young and very experimental app coming with absolutely no warranties." type="warning" show-icon)
  .form

    h4.formlabel Choose your input

    el-select(v-model='form.input_type')
      el-option(value='overhangs' label='A list of overhangs')
      el-option(value='constructs' label='Final constructs generated with CUBA')
      el-option(value='parts'  label='Individual parts (and backbone)')
    
    div(v-if="form.input_type === 'overhangs'")
      h4.formlabel Provide a list of overhangs
      .example(@click="form.overhangs_list = 'ATTG, GCTA, ACAT, CAAA, CCAG'"
             style='cursor: pointer; color: grey; margin: 0.5em;') Example
      el-input(type='textarea', :rows='4'
               v-model='form.overhangs_list',
               placeholder='Enter comma-separated overhangs, e.g. "ATTG, TTAC, ..."')
    
    div(v-if="form.input_type === 'constructs'")
      
      h4.formlabel Provide assembly genbank records

      collapsible(title='Examples')
        file-example(filename='example_assemblies.zip',
                     fileHref='/static/file_examples/predict_bad_clone_rates/example_assemblies.zip',
                     @input="function (e) {form.constructs_files.push(e); form.backbone_string='Resistance'}",
                     imgSrc='/static/file_examples/generic_logos/sequences_records.png')
          p.
            Genbank of constructs from the EMMA standard, generated using the CUBA
            "simulate Golden Gate" app (this gives them annotation which allow to
            detect overhangs).
      //- Accepted formats: FASTA, Genbank, Snapgene
      files-uploader(v-model='form.constructs_files', :multiple='true',
                     tip="Accepted formats: FASTA, Genbank, Snapgene")
      p Backbone indicator:
        el-input(v-model='form.backbone_string', placeholder="Enter a text present in the backbone's annotations")

    div(v-if="form.input_type === 'parts'")

      h4.formlabel Provide parts to assemble

        helper Upload parts sequences for one or several assemblies. Don't forget the receptor vector(s).
      collapsible(title='Examples')
        file-example(filename='example_genetic_parts_and_backbone.zip',
                     fileHref='/static/file_examples/simulate_gg_assemblies/example_genetic_parts_and_backbone.zip',
                     @input='function (e) {form.parts.push(e)}',
                     imgSrc='/static/file_examples/simulate_gg_assemblies/example_genetic_parts.png')
          p.
            Genbank records of five parts (A, A2, B, B2, C) and receptor vector. the parts can go into
            one of 3 possible slots, forming a total of four possible assemblies.
        file-example(filename='parts_missing_connectors.zip',
                    fileHref='/static/file_examples/simulate_gg_assemblies/parts_missing_connectors.zip',
                    @input='function (e) {form.parts.push(e)}',
                    imgSrc='/static/file_examples/simulate_gg_assemblies/example_genetic_parts.png')
          p.
            Genbank records of 11 parts which could form a circular assembly if they were completed
            with connectors. Tick "autoselect connectors" below and load the "example connectors",
            then press "Predict final constructs".
      //- Accepted formats: FASTA, Genbank, Snapgene
      files-uploader(v-model='form.parts_files', :multiple='true',
                     tip="Accepted formats: FASTA, Genbank, Snapgene")
      
      p: el-checkbox(v-model='form.use_assembly_plan') Provide a list of assemblies
      .use-assembly-plan(v-if='form.use_assembly_plan').animated.flipInX

        h4.formlabel Provide an assembly list

          helper.
            Only the connector parts necessary to obtain assemblies will be
            selected and added to the other parts.
        collapsible(title='Examples')
          file-example(filename='example_assembly_plan.xls',
                      fileHref='/static/file_examples/simulate_gg_assemblies/example_assembly_plan.xls',
                      @input='function (e) {form.assembly_plan = e}',
                      imgSrc='/static/file_examples/generic_logos/spreadsheet.svg')
            p.
              A picklist using the "example genetic parts" given as examples above.
        files-uploader(v-model='form.assembly_plan', :multiple='false')
        p Part names in the picklist refer to <br/>
          el-select(v-model='form.use_file_names_as_ids')
            el-option(:value='false', label="The IDs of the provided records")
            el-option(:value='true', label="The file names (without extension)")
        
    h4.formlabel Dataset parameters

    p Use data at Temp./Duration:
      el-select(v-model='form.dataset_tm', style='width: 100px; margin-left: 10px', size='small')
        el-option(value='25C' label='25C')
        el-option(value='37C' label='37C')
      el-select(v-model='form.dataset_duration', style='width: 100px;', size='small')
        el-option(value='01h' label='1h')
        el-option(value='18h' label='18h')
    p Dataset correction parameter
      el-input-number(v-model='form.dataset_corrective_factor', :min='0',
                      :max='3', :step='0.05', size='small',
                      style='width: 150px; ; margin-left: 10px')
    p In plots, mask interactions with value under
      el-input-number(v-model='form.rate_limit_in_plots', :min='0',
                      :max='10000', :step='1', size='small',
                      style='width: 150px; ; margin-left: 10px')

    backend-querier(:form='form', :backendUrl='infos.backendUrl',
                    :validateForm='validateForm', submitButtonText='Analyze',
                    v-model='queryStatus')
    el-alert(v-if='!queryStatus.polling.inProgress && queryStatus.requestError',
             :title="queryStatus.requestError", type="error", :closable="false")
  .results(v-if='!queryStatus.polling.inProgress && queryStatus.result.pdf_data')
    iframe(:src="queryStatus.result.pdf_data" 
            style="width: 90%; max-width: 1000px; height:800px;"
            frameborder="0")
    download-button(v-if='queryStatus.result.csv_file_data',
                    text='Download Summary Spreadsheet'
                    :filedata='queryStatus.result.csv_file_data')
  powered-by(:softwareNames='infos.poweredby')
</template>

<script>
import learnmore from '../../components/widgets/LearnMore'
import sequencesuploader from '../../components/widgets/SequencesUploader'
import digestionset from '../../components/widgets/DigestionSelectorSet'
import ladderselector from '../../components/widgets/LadderSelector'

var infos = {
  title: 'Predict Bad Clone Rates',
  navbarTitle: 'Predict Bad Clone Rates',
  path: 'predict-bad-clone-rates',
  description: '',
  backendUrl: 'start/predict_bad_clone_rates',
  icon: require('../../assets/images/predict_bad_clone_rates.svg'),
  poweredby: ['kappagate']
}

export default {
  data () {
    return {
      form: {
        input_type: 'overhangs',
        constructs_files: [],
        overhangs_list: '',
        parts_files: [],
        use_assembly_plan: false,
        assembly_plan: null,
        use_file_names_as_ids: false,
        rate_limit_in_plots: 100,
        dataset_tm: '25C',
        dataset_duration: '18h',
        dataset_corrective_factor: 1.0,
        backbone_string: ''
      },
      infos: infos,
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
    digestionset,
    ladderselector
  },
  infos: infos,
  methods: {
    handleSuccess (evt) {
      console.log(evt)
    },
    validateForm () {
      var errors = []
      if ((this.form.input_type === 'overhangs') && (this.form.overhangs_list.split(',').length < 3)) {
        errors.push('Provide at least 3 overhangs')
      }
      if ((this.form.input_type === 'parts') && (this.form.parts_files.length < 3)) {
        errors.push('Provide at least 3 parts')
      }
      if ((this.form.input_type === 'constructs') && (this.form.constructs_files.length === 0)) {
        errors.push('Provide some construct files')
      }
      return errors
    }
  }
}
</script>

<style lang='css' scoped>

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
</style>
