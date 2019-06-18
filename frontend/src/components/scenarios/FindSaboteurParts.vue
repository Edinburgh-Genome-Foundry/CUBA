<template lang="pug">
.page
  h1  {{ infos.title }}
  img.icon.center-block(slot='title-img', :src='infos.icon')
  p.scenario-description.
    Find which parts are associated with assembly failure!
  web-links(:emailSubject="'[CUBA] Feedback on web app: ' + infos.title",
            tweetMessage="Sequenticons are human-friendly, visual DNA sequence identifiers.",
            :tweetUrl="'https://cuba.genomefoundry.org/' + infos.path")
  //- learnmore Bla bla bla

  .form

    h4.formlabel Method
    el-select(v-model='form.method')
      el-option(value='statistical' label='Use a statistical method')
      el-option(value='logical' label='Use a logical method')
      

    h4.formlabel Assemblies data for {{form.method}} method
    collapsible(title='Examples')
      file-example(:filename='`${form.method}_assembly_data.csv`',
                   :fileHref='`/static/file_examples/find_saboteur_parts/${form.method}_assembly_data.csv`',
                   @input='function (e) { form.assemblies_data_file = e }'
                   imgSrc='/static/file_examples/generic_logos/spreadsheet.svg')
    files-uploader(v-model='form.assemblies_data_file', text="Drop files (or click to select)",
                  help='CSV files only :)', :multiple='false')
    h4.formlabel Parts-per-construct data
    el-select(v-model='form.parts_input_type')
      el-option(value='in_data' label='It is provided in the above spreadsheet')
      el-option(value='assembly_plan' label="I'll provide it as an assembly plan")
      el-option(value='records' label="I'll provide final constructs records")
    div(v-if="form.parts_input_type !== 'in_data'")
      h4.formlabel(v-if="form.parts_input_type === 'assembly_plan'") Assembly plan
      h4.formlabel(v-if="form.parts_input_type === 'records'") Assembly records
      files-uploader(v-model='form.parts_data_files', text="Drop files (or click to select)",
                    help='CSV files only :)', :multiple='true')
    backend-querier(:form='form', :backendUrl='infos.backendUrl',
                    :validateForm='validateForm', submitButtonText='Find saboteurs',
                    v-model='queryStatus')
    el-alert(v-show='queryStatus.requestError', :title="queryStatus.requestError",
       type="error", :closable="false")

  .results(v-if='!queryStatus.polling.inProgress && Object.keys(queryStatus.result).length')
    download-button(v-if='queryStatus.result.report',
                    text='Download the report',
                    :filedata='queryStatus.result.report')
    .textual-report(v-if='queryStatus.result.suspicious')
      div(v-if='queryStatus.result.saboteurs.length')
        h3 Bad parts identified !
        p.
          The following parts appear only in failed assemblies and for each of
          them there is at least one assembly in which all other parts
          appear in successful assemblies. Therefore these parts are
          #[b probably bad].
          :
        ul
          li(v-for='partName in queryStatus.result.suspicious', :key='partName') {{partName}}
      div(v-else)
        p No saboteur parts were found
      div(v-if='queryStatus.result.suspicious.length')
        h3 Suspicious parts identified !
        p.
          The following parts appear only in failed assemblies but are #[b only suspicious]
          as there is no other evidence that they are bad (some of them may be,
          some may not be).
        ul
          li(v-for='partName in queryStatus.result.suspicious', :key='partName') {{partName}}


  powered-by(:softwareNames='infos.poweredby')
</template>

<script>
import learnmore from '../../components/widgets/LearnMore'

var infos = {
  title: 'Find Saboteur Parts',
  navbarTitle: 'Find Saboteur Parts',
  path: 'find_saboteur_parts',
  description: '',
  backendUrl: 'start/find_saboteur_parts',
  icon: require('../../assets/images/find_saboteur_parts.svg'),
  poweredby: ['saboteurs']
}

export default {
  data () {
    return {
      infos: infos,
      form: {
        method: 'statistical',
        assemblies_data_file: null,
        parts_input_type: 'in_data',
        parts_data_files: []
      },
      queryStatus: {
        polling: {},
        result: {},
        requestError: ''
      }
    }
  },
  components: {
    learnmore
  },
  infos: infos,
  methods: {
    validateForm () {
      var errors = []
      if (!this.form.assemblies_data_file) {
        errors.push('Provide files !')
      }
      return errors
    }
  }
}
</script>

<style lang='scss' scoped>

.form {
  margin: 50px auto;
  max-width: 500px;
}

h4.formlabel {
  text-align: center;
  text-transform: uppercase;
  margin-top: 40px
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
