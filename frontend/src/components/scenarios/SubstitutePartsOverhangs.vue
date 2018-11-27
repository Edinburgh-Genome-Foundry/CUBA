<template lang="pug">

.page
  h1 {{infos.title}}
  web-links(:emailSubject="'[CUBA] Feedback on web app: ' + infos.title",
            tweetMessage="Substitute overhangs in parts",
            :tweetUrl="'https://cuba.genomefoundry.org/' + infos.path")
  img.icon.center-block(slot='title-img', :src='infos.icon')
  p.scenario-description.
    Use this app to easily change the position (a.k.a "slot") of parts in
    Golden-Gate like assembly standards.


  .form
    h4.formlabel Provide parts to modify
    collapsible(title='Examples')
      file-example(filename='example_genetic_parts_and_backbone.zip',
                   fileHref='/static/file_examples/simulate_gg_assemblies/example_genetic_parts_and_backbone.zip',
                   @input='function (e) {form.parts.push(e)}',
                   imgSrc='/static/file_examples/simulate_gg_assemblies/example_genetic_parts.png')
        p.
          Genbank records of five parts (A, A2, B, B2, C) and receptor vector.
    files-uploader(v-model='form.parts', :multiple='true',
                   tip="Accepted formats: FASTA, Genbank, Snapgene")
    el-checkbox(v-model='form.use_file_names_as_ids') Use file names as IDs.

    h4.formlabel Select an enzyme
    .enzymes-radio
      el-row(:gutter='60')
        el-col(v-for='enzyme in enzymes', :md='6', :sm='6', :xs='24', :key='enzyme')
          el-radio(v-model='form.enzyme', class="radio", :label='enzyme') {{enzyme}}

    h4.formlabel Substitutions
    el-input(type='textarea' :autosize='true', v-model='form.substitutions',
             placeholder='One substitution per line, e.g. ATGC => CTCA')
    p
      el-checkbox(v-model='form.return_linear_parts') Return a FASTA of linear parts for ordering.

    backend-querier(:form='form', :backendUrl='infos.backendUrl',
                    :validateForm='validateForm', submitButtonText='Substitute overhangs',
                    v-model='queryStatus')
    el-alert(v-if='queryStatus.requestError', :title="queryStatus.requestError",
       type="error", :closable="false")
  .results(v-if='!queryStatus.polling.inProgress  && queryStatus.result.file')
    download-button(v-if='queryStatus.result.file',
                    :filedata='queryStatus.result.file')
    center
      iframe(:src="queryStatus.result.pdf_report.data" 
              style="width: 90%; max-width: 800px; height:500px;" frameborder="0")
      //- .unmodified(v-if='queryStatus.result.unmodified.length')
      //-   p These parts were unmodified: {{ queryStatus.result.unmodified.join(' ') }}
  powered-by(:softwareNames='infos.poweredby')
</template>

<script>
import learnmore from '../../components/widgets/LearnMore'
import sequencesuploader from '../../components/widgets/SequencesUploader'
import digestionset from '../../components/widgets/DigestionSelectorSet'
import ladderselector from '../../components/widgets/LadderSelector'

var infos = {
  title: 'Substitute Parts Overhangs',
  navbarTitle: 'Substitute Parts Overhangs',
  path: 'substitute_parts_overhangs',
  description: '',
  backendUrl: 'start/substitute_parts_overhangs',
  icon: require('../../assets/images/substitute_parts_overhangs.svg'),
  poweredby: ['dnacauldron']
}

export default {
  data () {
    return {
      form: {
        parts: [],
        substitutions: '',
        enzyme: 'Autoselect',
        use_file_names_as_ids: true,
        return_linear_parts: true
      },
      enzymes: ['BsaI', 'BsmBI', 'BbsI', 'Autoselect'],
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
      if (this.form.parts.length === 0) {
        errors.push('Provide parts')
      }
      if (this.form.substitutions.length === 0) {
        errors.push('Provide some substitutions')
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
