<template lang="pug">
div
  h1  {{ infos.title }}
  img.icon.center-block(slot='title-img', :src='infos.icon')
  p.center.
    Submit parts and a receptor vector. Get an annotated Genbank of the
    resulting construct(s).
  web-links(:emailSubject="'[CUBA] Feedback on web app: ' + infos.title",
            tweetMessage="Simple online cloning simulator for Golden Gate:",
            :tweetUrl="'http://cuba.genomefoundry.org/' + infos.path")

  .form
    h4.formlabel Select an enzyme
    .enzymes-radio
      el-row(:gutter='60')
        el-col(v-for='enzyme in enzymes', :md='6', :sm='6', :xs='24', :key='enzyme')
          el-radio(v-model='form.enzyme', class="radio", :label='enzyme') {{enzyme}}

    h4.formlabel Provide Parts and a receptor vector
    sequencesuploader(v-model='form.parts', :multiple='true',
                      text="Drop multiple Genbank/Fasta (or click to select)")
    h4.formlabel (Optional) Provide connectors
    .select-connectors
        p.
          Only the connector parts necessary to obtain assemblies will be
          selected and added to the other parts.
        sequencesuploader(v-model='form.connectors', :multiple='true',
                          text="Drop multiple Genbank/Fasta (or click to select)")


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
  backendUrl: 'start/simulate_cloning',
  icon: require('assets/images/assembly.svg'),
  poweredby: ['dnacauldron', 'dnafeaturesviewer']
}

export default {
  data: function () {
    return {
      enzymes: ['BsaI', 'BsmBI', 'BbsI', 'Autoselect'],
      form: {
        enzyme: 'Autoselect',
        parts: [],
        connectors: []
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
    handleSuccess: function (evt) {
      console.log(evt)
    },
    validateForm: function () {
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
