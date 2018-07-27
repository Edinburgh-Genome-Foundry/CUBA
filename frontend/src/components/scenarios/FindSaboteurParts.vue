<template lang="pug">
.page
  h1  {{ infos.title }}
  img.icon.center-block(slot='title-img', :src='infos.icon')
  p.scenario-description Generate the #[a(href='https://github.com/Edinburgh-Genome-Foundry/sequenticon') sequenticons] corresponding to your sequences.
  web-links(:emailSubject="'[CUBA] Feedback on web app: ' + infos.title",
            tweetMessage="Sequenticons are human-friendly, visual DNA sequence identifiers.",
            :tweetUrl="'http://cuba.genomefoundry.org/' + infos.path")
  //- learnmore Bla bla bla

  .form

    h4.formlabel Assemblies data
    collapsible(title='Examples')
      file-example(filename='assemblies_data.csv',
                   fileHref='/static/file_examples/find_saboteur_parts/assemblies_data.csv',
                   @input='function (e) { form.assemblies_data_file = e }'
                   imgSrc='/static/file_examples/generic_logos/spreadsheet.svg')
    files-uploader(v-model='form.assemblies_data_file', text="Drop files (or click to select)",
                  help='Fasta or genbank. No file too large please :)', :multiple='false')

    backend-querier(:form='form', :backendUrl='infos.backendUrl',
                    :validateForm='validateForm', submitButtonText='Evaluate',
                    v-model='queryStatus')
    el-alert(v-show='queryStatus.requestError', :title="queryStatus.requestError",
       type="error", :closable="false")

  .results(v-if='!queryStatus.polling.inProgress && queryStatus.result')
    download-button(v-if='queryStatus.result.report',
                    text='Download the report',
                    :filedata='queryStatus.result.report')
    .images(v-if="queryStatus.result.images")
      el-card.box-card(v-for="name_tag in queryStatus.result.images"
                       :key='name_tag[0]'
                       style="width:150px; display: inline-block; text-align: center; margin: 0.5cm;")
        .clearfix(slot="header")
          span {{ name_tag[0] }}
        .img(v-html="name_tag[1]")


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
        assemblies_data_file: null
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
