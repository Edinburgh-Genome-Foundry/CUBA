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
                   imgSrc='/static/file_examples/simulate_gg_assemblies/example_genetic_parts.png')
        p.
          Genbank records of five parts (A, A2, B, B2, C) and receptor vector. the parts can go into
          one of 3 possible slots, forming a total of four possible assemblies.
    sequencesuploader(v-model='form.parts', :multiple='true',
                      text="Drop multiple Genbank/Fasta (or click to select)")
    p: el-checkbox(v-model='form.use_assembly_list') Provide a list of assemblies
    .select-connectors(v-if='form.use_assembly_list').animated.flipInX
      h4.formlabel Provide an assembly list
        helper.
          Only the connector parts necessary to obtain assemblies will be
          selected and added to the other parts.")
      files-uploader(v-model='form.assembly_list')
    p: el-checkbox(v-model='form.select_connectors') Autoselect connectors  <br/>
    .select-connectors(v-if='form.select_connectors').animated.flipInX
      h4.formlabel Provide connectors
        helper.
          Connectors are neutral parts used to fill gaps in some assemblies.
          Drop here all the connector sequences you have, only the connectors
          necessary to obtain assemblies will be selected and added to the other parts.

      sequencesuploader(v-model='form.connectors', :multiple='true',
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
      .results-summary(v-if='queryStatus.result.preview',
                       v-html="queryStatus.result.preview.html")
  powered-by(:softwareNames='infos.poweredby')
</template>

<script>
import learnmore from '../../components/widgets/LearnMore'
import sequencesuploader from '../../components/widgets/SequencesUploader'
import digestionset from '../../components/widgets/DigestionSelectorSet'

var infos = {
  title: 'Simulate Golden Gate assemblies',
  navbarTitle: 'Simulate Golden Gate',
  path: 'simulate_gg_assemblies',
  description: '',
  backendUrl: 'start/simulate_gg_assemblies',
  icon: require('assets/images/assembly.svg'),
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
        use_assembly_list: false,
        assembly_list: null
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
