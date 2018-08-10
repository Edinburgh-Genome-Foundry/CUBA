<template lang="pug">
.page
  h1  {{ infos.title }}
  web-links(:emailSubject="'[CUBA] Feedback on web app: ' + infos.title",
            tweetMessage="Simple online cloning simulator for Golden Gate:",
            :tweetUrl="'http://cuba.genomefoundry.org/' + infos.path")
  img.icon.center-block(slot='title-img', :src='infos.icon')
  p.scenario-description.
    Submit parts and a receptor vector. Get an annotated Genbank of the
    resulting construct(s).


  .form
    h4.formlabel Select an enzyme
    .enzymes-radio
      el-row(:gutter='60')
        el-col(v-for='enzyme in enzymes', :md='6', :sm='6', :xs='24', :key='enzyme')
          el-radio(v-model='form.enzyme', class="radio", :label='enzyme') {{enzyme}}

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
    files-uploader(v-model='form.parts', :multiple='true',
                   tip="Accepted formats: FASTA, Genbank, Snapgene")
    //- sequencesuploader(v-model='form.parts', :multiple='true',
    //-                   text="Drop multiple Genbank/Fasta (or click to select)")
    p: el-checkbox(v-model='form.use_assembly_plan') Provide a list of assemblies
    .use-assembly-plan(v-if='form.use_assembly_plan').animated.flipInX
      h4.formlabel Provide an assembly list
        helper.
          Only the connector parts necessary to obtain assemblies will be
          selected and added to the other parts.")
      collapsible(title='Examples')
        file-example(filename='example_assembly_plan.xls',
                     fileHref='/static/file_examples/simulate_gg_assemblies/example_assembly_plan.xls',
                     @input='function (e) {form.assembly_plan = e}',
                     imgSrc='/static/file_examples/generic_logos/spreadsheet.svg')
          p.
            A picklist using the "example genetic parts" given as examples above.
      files-uploader(v-model='form.assembly_plan', :multiple='false')
      p: el-checkbox(v-model='form.single_assemblies') Ensure each line gives a single assembly
    p: el-checkbox(v-model='form.select_connectors') Autoselect connectors  <br/>
    .select-connectors(v-if='form.select_connectors').animated.flipInX
      h4.formlabel Provide connectors
        helper.
          Connectors are neutral parts used to fill gaps in some assemblies.
          Drop here all the connector sequences you have, only the connectors
          necessary to obtain assemblies will be selected and added to the other parts.
      collapsible(title='Examples')
        file-example(filename='emma_connectors.zip',
                     fileHref='/static/file_examples/simulate_gg_assemblies/emma_connectors.zip',
                     @input='function (e) {form.connectors.push(e)}',
                     imgSrc='/static/file_examples/generic_logos/spreadsheet.svg')
          p.
            Connectors part, some of which can complete the parts given above as an example of
            "parts missing connectors".

      files-uploader(v-model='form.connectors', :multiple='true',
                     text="Drop multiple Genbank/Fasta (or click to select)")
    p: el-checkbox(v-model='form.include_fragments') Include parts and fragments in report (slower)


    backend-querier(:form='form',
                    :backendUrl='infos.backendUrl',
                    :validateForm='validateForm',
                    submitButtonText='Predict final construct(s)',
                    v-model='queryStatus')
    el-alert(v-if='queryStatus.requestError && !queryStatus.polling.inProgress',
             :title="queryStatus.requestError",
             type="error",
             :closable="false")
    .results(v-if='!queryStatus.polling.inProgress')
      download-button(v-if='queryStatus.result.file',
                      :filedata='queryStatus.result.file')
      .div(v-if='queryStatus.result.infos')
        .nConstructs(v-if='queryStatus.result.infos.nconstructs || queryStatus.result.infos.nconstructs === 0')
          span(v-if='queryStatus.result.infos.nconstructs === 0' style='color: red').
            No construct could be found. See report for more.
          span(v-else) {{queryStatus.result.infos.nconstructs}} constructs generated !

        .errors(v-if='queryStatus.result.infos.errors && queryStatus.result.infos.errors.length')
          p The following errors occured in the assembly plan (see report for more)
          ul
            li(v-for='[construct, error] in queryStatus.result.infos.errors') #[b {{construct}}]: {{error}}

      .didyoumean(v-if='queryStatus.result.unknown_parts')
        p There were parts in the assembly plan without any correspondance in the records:
        .div(v-for='didyoumean, name in queryStatus.result.unknown_parts')
          p Part {{name}}: did you mean...
          ul
            li(v-for='candidate in didyoumean') {{candidate}}

      .results-summary(v-if='queryStatus.result.preview',
                       v-html="queryStatus.result.preview.html")
  powered-by(:softwareNames='infos.poweredby')
</template>

<script>
import learnmore from '../../components/widgets/LearnMore'
import sequencesuploader from '../../components/widgets/SequencesUploader'
import digestionset from '../../components/widgets/DigestionSelectorSet'

var infos = {
  title: 'Simulate Golden Gate Assemblies',
  navbarTitle: 'Simulate Golden Gate',
  path: 'simulate_gg_assemblies',
  description: '',
  backendUrl: 'start/simulate_gg_assemblies',
  icon: require('../../assets/images/assembly.svg'),
  poweredby: ['dnacauldron', 'dnafeaturesviewer']
}

export default {
  data () {
    return {
      enzymes: ['BsaI', 'BsmBI', 'BbsI', 'Autoselect'],
      form: {
        enzyme: 'Autoselect',
        parts: [],
        connectors: [],
        select_connectors: false,
        include_fragments: false,
        use_assembly_plan: false,
        assembly_plan: null,
        single_assemblies: true
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
    sequencesuploader,
    learnmore,
    digestionset
  },
  infos: infos,
  methods: {
    handleSuccess (evt) {
      console.log(evt)
    },
    validateForm () {
      var errors = []
      if (this.form.parts.length === 0) {
        errors.push('Provide at least 2 files: one or more parts and a receptor.')
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

.enzymes-radio {
  width: 400px;
  margin: 0 auto;
  .el-radio {
    margin-bottom: 20px;
  }
}
</style>
