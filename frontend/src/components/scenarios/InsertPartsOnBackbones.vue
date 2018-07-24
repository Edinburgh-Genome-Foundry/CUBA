<template lang="pug">
.page
  h1  {{ infos.title }}
  web-links(:emailSubject="'[CUBA] Feedback on web app: ' + infos.title",
            tweetMessage="Automatically put part sequences on a backbone !",
            :tweetUrl="'http://cuba.genomefoundry.org/' + infos.path")
  img.icon.center-block(slot='title-img', :src='infos.icon')
  p.scenario-description.
    Quickly swap a part (from a golden-gate-like based assembly standard) from
    one vector to another (with same overhangs).
    Provide the donor vector to insert in, and a list of insert parts
    (either as linearized or plasmid records). Get a version of your part(s)
    inserted in the donor vector.


  .form
    h4.formlabel You want to

    el-select(v-model='form.mode', placeholder='Select')
      el-option(label='Put one or several parts on one backbone', value='one_backbone')
      el-option(label='Auto-select backbones for different parts', value='select_backbone')


    h4.formlabel Select an enzyme
    .enzymes-radio
      el-row(:gutter='60')
        el-col(v-for='enzyme in enzymes', :md='6', :sm='6', :xs='24', :key='enzyme')
          el-radio(v-model='form.enzyme', class="radio", :label='enzyme') {{enzyme}}

    div(v-if="form.mode === 'one_backbone'")
      h4.formlabel Provide a backbone
      files-uploader(v-model='form.backbone', :multiple='false',
                     tip="Accepted formats: Genbank/Fasta/txt/.doc")
      h4.formlabel Provide parts to insert
          helper The parts can be either on a linear or plasmid sequence.

    div(v-if="form.mode === 'select_backbone'")
      h4.formlabel Provide candidate backbones
      files-uploader(v-model='form.backbones', :multiple='true',
                     tip="Accepted formats: Genbank/Fasta/txt/.doc")
      h4.formlabel Provide parts to insert
          helper The parts can be either on a linear or plasmid sequence.


    files-uploader(v-model='form.inserts', :multiple='true',
                   tip="Accepted formats: Genbank/Fasta/txt/.doc")


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
  title: 'Insert parts on backbones',
  navbarTitle: 'Insert parts on backbones',
  path: 'insert_parts_on_backbones',
  description: '',
  backendUrl: 'start/insert_parts_on_backbones',
  icon: require('../../assets/images/insert_parts_on_backbones.svg'),
  poweredby: ['dnacauldron', 'dnafeaturesviewer']
}

export default {
  data () {
    return {
      enzymes: ['BsaI', 'BsmBI', 'BbsI', 'Autoselect'],
      form: {
        enzyme: 'Autoselect',
        backbone: null,
        backbones: [],
        inserts: [],
        mode: 'one_backbone'
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
      if (this.form.inserts.length === 0) {
        errors.push('Provide at least 1 insert.')
      }
      if (this.form.donor_vector === {}) {
        errors.push('Provide at least 1 insert.')
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

.form, .page p {
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
