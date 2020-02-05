<template lang="pug">
.page
  h1  {{ infos.title }}
  web-links(:emailSubject="'[CUBA] Feedback on web app: ' + infos.title",
            tweetMessage="Simple online cloning simulator for Golden Gate:",
            :tweetUrl="'https://cuba.genomefoundry.org/' + infos.path")
  img.icon.center-block(slot='title-img', :src='infos.icon')
  p.scenario-description.
    Submit parts and a receptor vector. Get an annotated Genbank of the
    resulting construct(s). #[b Note:] if you want 


  .form

    h4.formlabel Provide an assembly plan
      helper.
        Only the connector parts necessary to obtain assemblies will be
        selected and added to the other parts.
    collapsible(title='Examples')
      file-example(filename='multi_assembly.xlsx',
                    fileHref='/static/file_examples/simulate_multi_method_assemblies/multi_assembly.xlsx',
                    @input='function (e) {form.assembly_plan = e}',
                    imgSrc='/static/file_examples/generic_logos/spreadsheet.svg')
        p.
          An assembly plan with many different cloning methods used (Golden Gate,
          Gibson, BioBrick, LCR and Basic assembly).
    files-uploader(v-model='form.assembly_plan', :multiple='false')
    p Part names in the picklist refer to <br/>
      el-select(v-model='form.use_file_names_as_ids')
        el-option(:value='false', label="The IDs of the provided records")
        el-option(:value='true', label="The file names (without extension)")


    h4.formlabel Provide parts to assemble
      helper Upload all the parts sequences for your assemblies. Don't forget the receptor vectors.
    collapsible(title='Examples')
      file-example(filename='part_sequences.zip',
                   fileHref='/static/file_examples/simulate_multi_method_assemblies/part_sequences.zip',
                   @input='function (e) {form.parts.push(e)}',
                   imgSrc='/static/file_examples/simulate_gg_assemblies/example_genetic_parts.png')
        p.
          Genbank and fasta records of all parts needed in the assembly plan
    //- Accepted formats: FASTA, Genbank, Snapgene
    files-uploader(v-model='form.parts', :multiple='true',
                   tip="Accepted formats: FASTA, Genbank, Snapgene")
    //- sequencesuploader(v-model='form.parts', :multiple='true',
    //-                   text="Drop multiple Genbank/Fasta (or click to select)")

    p
      el-select(v-model='form.topology' size='small')
        el-option(value='circular' label='All sequences are circular')
        el-option(value='linear' label='All sequences are linear')
        el-option(value='default_to_circular' label='Autodetect each sequence\'s topology  (default to circular)')
        el-option(value='default_to_linear' label='Autodetect each sequence\'s topology (default to linear)')
    

    p: el-checkbox(v-model='form.select_connectors') Autoselect connectors  <br/>
    .select-connectors(v-if='form.select_connectors').animated.flipInX
      h4.formlabel Provide connectors
        helper.
          Connectors are neutral parts used to fill gaps in some assemblies.
          Drop here all the connector sequences you have, only the connectors
          necessary to obtain assemblies will be selected and added to the other parts.

      p.
        For each "connector_collection" referenced in your assembly plan, for
        instance "CC1" in the example assembly plan, gather all the connectors in
        the collection under a same-named zip file or fasta file
        (CC1.zip or CCA.file).

      collapsible(title='Examples')
        file-example(filename='CC1.zip',
                     fileHref='/static/file_examples/simulate_multi_method_assemblies/CC1.zip',
                     @input='function (e) {form.connectors.push(e)}',
                     imgSrc='/static/file_examples/generic_logos/sequence_file.svg')
          p.
            Connectors part, some of which can complete assembly one assembly
            in the example assembly plan.

      files-uploader(v-model='form.connectors', :multiple='true',
                     text="Drop multiple Genbank/Fasta (or click to select)")

    p
      el-select(v-model='form.include_fragment_plots')
        el-option(value='on_failure' label="Plot parts and fragments on assembly failure only")
        el-option(value='yes' label="Plot parts and fragments (slower)")
        el-option(value='no' label="Don't plot parts and fragments")
    p
      el-select(v-model='form.include_graph_plots')
        el-option(value='on_failure' label="Plot graphs on assembly failure only")
        el-option(value='yes' label="Always plot graph (slower)")
        el-option(value='no' label="Don't plot graph")


    backend-querier(:form='form',
                    :backendUrl='infos.backendUrl',
                    :validateForm='validateForm',
                    submitButtonText='Predict final construct(s)',
                    v-model='queryStatus')
    progress-bars(:bars='queryStatus.polling.data.bars', :order="['assembly']"
                  v-if='queryStatus.polling.inProgress && queryStatus.polling.data')
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
      .stats(v-if='queryStatus.result.assembly_stats')
        ul
          li #[b Valid assemblies:] {{queryStatus.result.assembly_stats.valid_assemblies}}
          li #[b Errored assemblies:] {{queryStatus.result.assembly_stats.errored_assemblies}}
          li #[b Cancelled assemblies:] {{queryStatus.result.assembly_stats.cancelled_assemblies}}
      .stats(v-else) {{queryStatus.result.n_constructs}} constructs generated!
      .errors(v-if='queryStatus.result.errors && queryStatus.result.errors.length')
        p The following errors occured in the assembly plan (see report for more)
        ul
          li(v-for='error in queryStatus.result.errors') {{error}}
      
  powered-by(:softwareNames='infos.poweredby')
</template>

<script>
import learnmore from '../../components/widgets/LearnMore'
import sequencesuploader from '../../components/widgets/SequencesUploader'
import digestionset from '../../components/widgets/DigestionSelectorSet'

var infos = {
  title: 'Simulate multi-method assemblies',
  navbarTitle: 'Simulate multi-method assemblies',
  path: 'simulate_multi_method_assemblies',
  description: '',
  backendUrl: 'start/simulate_multi_method_assemblies',
  icon: require('../../assets/images/simulate_multi_method_assembly.svg'),
  poweredby: ['dnacauldron', 'dnafeaturesviewer']
}

export default {
  data () {
    return {
      form: {
        enzyme: 'auto',
        parts: [],
        connectors: [],
        select_connectors: false,
        include_fragment_plots: 'on_failure',
        include_assembly_plots: true,
        include_graph_plots: 'on_failure',
        assembly_plan: null,
        use_file_names_as_ids: false,
        topology: 'default_to_circular'
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
        errors.push('Provide at least some parts.')
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
