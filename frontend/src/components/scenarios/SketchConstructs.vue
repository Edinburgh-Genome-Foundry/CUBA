<template lang="pug">
.page
  h1  {{ infos.title }}
  web-links(:emailSubject="'[CUBA] Feedback on web app: ' + infos.title",
            tweetMessage="Simple online DNA construct sketcher:",
            :tweetUrl="'https://cuba.genomefoundry.org/' + infos.path")
  img.icon.center-block(slot='title-img', :src='infos.icon')
  p.scenario-description.
    Sketch one or several assemblies below, then export as PDF or PNG.

  .sketcher
    sketcher(v-model='form.sketchesData')
  .form
    el-form(:model='form' label-width="120px")
      el-form-item(label='Format')
        el-select(v-model='form.format')
          el-option(value='PDF' label='PDF')
          el-option(value='jpeg' label='JPEG')
          el-option(value='PNG' label='PNG')
        //- el-radio(v-model='form.format' label='jpeg') JPEG
        //- el-radio(v-model='form.format' label='PNG') PNG
      el-form-item(label='Orientation' v-if="form.format === 'PDF'")
        el-select(v-model='form.orientation')
             el-option(value='portrait' label='Protrait')
             el-option(value='landscape' label='Landscape')
      el-form-item(label='Width' v-else)
        el-input-number(v-model='form.width', :min='50', :max='1000' size='small')
      el-form-item(label='Font')
        el-select(v-model='form.font')
          el-option.open-sans(value='Open Sans' label='Open Sans')
          el-option.open-sans(value='Raleway' label='Raleway')
          el-option.open-sans(value='Cabin sketch' label='Cabin sketch')
      el-form-item(label='Font size')
        el-input-number(v-model='form.fontsize', :min='5', :max='18' size='small')




  backend-querier(:form='form', :backendUrl='infos.backendUrl',
                  :validateForm='validateForm', submitButtonText='Render',
                  v-model='queryStatus')
  el-alert(v-if='queryStatus.requestError && !queryStatus.polling.inProgress',
           :title="queryStatus.requestError",
     type="error", :closable="false")
  .results(v-if='!queryStatus.polling.inProgress')
    download-button(v-if='queryStatus.result.file',
                    text='Download',
                    :filedata='queryStatus.result.file')
    .results-summary(v-if='queryStatus.result.preview',
                     v-html="queryStatus.result.preview.html")
  powered-by(:softwareNames='infos.poweredby')
</template>

<script>
import sketcher from '../widgets/SeqSketcher/SeqSketcher.vue'
import tools from '../../tools'

var infos = {
  title: 'Sketch Constructs',
  navbarTitle: 'Sketch Constructs',
  path: 'sketch_constructs',
  description: '',
  backendUrl: 'start/sketch_constructs',
  icon: require('../../assets/images/sketch_constructs.svg'),
  poweredby: ['caravagene']
}

export default {
  data () {
    return {
      form: {
        format: 'PDF',
        orientation: 'portrait',
        fontsize: 14,
        width: 600,
        font: 'Raleway',
        sketchesData: {
          title: '',
          note: '',
          constructs: [
            {
              name: '',
              note: '',
              id: tools.generateRandomID(),
              parts: [
                {
                  category: 'promoter',
                  label: 'prom1',
                  sublabel: 'Hover us !',
                  id: 0
                },
                {
                  category: 'CDS',
                  label: 'gene with a very very long name',
                  id: 1
                },
                {
                  category: 'terminator',
                  label: 'Term1',
                  subscript: 'Some subscript',
                  id: 2
                }
              ]
            }
          ]
        }
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
    sketcher
  },
  infos: infos,
  methods: {
    handleSuccess (evt) {
      console.log(evt)
    },
    validateForm () {
      var errors = []
      return errors
    }
  }
}
</script>

<style lang='scss' scoped>

.page {
  width: 1500px;
}

h4.formlabel {
  text-align: center;
  text-transform: uppercase;
  margin-top: 40px
}

.form {
  .raleway {font-family: 'Raleway'}
  .open-sans {font-family: 'Open Sans'}
  .cabin-sketch {font-family: 'Cabin Sketch'}
  width: 500px;
  margin: 0 auto;
}
.sketcher {
  margin-top: 2em;
  padding: 0.5em;
  padding-top: 1em;
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
  width: 99%
}

.enzymes-radio {
  width: 400px;
  margin: 0 auto;
  .el-radio {
    margin-bottom: 20px;
  }
}
</style>
