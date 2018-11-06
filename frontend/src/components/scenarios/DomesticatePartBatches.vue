<template lang="pug">
.page
  h1 {{infos.title}}
  web-links(:emailSubject="'[CUBA] Feedback on web app: ' + infos.title",
            tweetMessage="Domesticate parts for EMMA and other standards:",
            :tweetUrl="'https://cuba.genomefoundry.org/' + infos.path")
  img.icon.center-block(slot='title-img', :src='infos.icon')
  p.scenario-description.
    Provide parts and have them domesticated for EMMA assembly or another standard.

  el-alert(title="This app is a stub. Not all features may work properly" type="warning" show-icon)

  .form
    //- collapsible(title='Examples')
    //-   file-example(filename='emma_parts.zip',
    //-                @input='function (e) {form.parts_infos.push(e)}',
    //-                fileHref='/static/file_examples/create_assembly_picklists/emma_parts.zip',
    //-                imgSrc='/static/file_examples/generic_logos/sequences_records.png')
    //-     p.
    //-       ZIP archive containing the sequences of standard EMMA parts.
    h4.formlabel Standard

    p
      el-select(v-model='form.standard')
        el-option(value='EMMA', label='EMMA')
        el-option(value='custom', label='Custom')
    el-form(v-if="form.standard == 'custom'")
      el-form-item(label='Name')
        el-input(v-model='form.standard_name', placeholder='Name of the standard')
      el-form-item(label='Definition')
        files-uploader(v-model='form.standard_definition', tip='csv or excel', :multiple='false')
    hr
    h4.formlabel Parts
    files-uploader(v-model='form.parts', tip='Genbank/Fasta/Snapgene sequences, or spreadsheet')

    backend-querier(:form='form', :backendUrl='infos.backendUrl',
                    :validateForm='validateForm', submitButtonText='Domesticate',
                    v-model='queryStatus')
    progress-bars(:bars='queryStatus.polling.data.bars', :order="['record', 'primer']"
                  v-if='queryStatus.polling.inProgress && queryStatus.polling.data')
    el-alert(v-if='queryStatus.requestError', :title="queryStatus.requestError",
       type="error", :closable="false")

    .results(v-if='!queryStatus.polling.inProgress')
      p(v-if='queryStatus.result.nfails > 0') There were {{queryStatus.result.nfails}} errors, see report for more.
      p(v-else) Everything seems to be fine. See report.
      download-button(v-if='queryStatus.result.file',
        text='Download Report',
        :filedata='queryStatus.result.file')
  powered-by(:softwareNames='infos.poweredby')
</template>

<script>

var infos = {
  title: 'Domesticate Part Batches',
  navbarTitle: 'Domesticate Part Batches',
  path: 'domesticate_part_batches',
  description: '',
  backendUrl: 'start/domesticate_part_batches',
  icon: require('../../assets/images/domesticate_part_batches.svg'),
  poweredby: ['genedom']
}

export default {
  data () {
    return {
      form: {
        parts: [],
        standard: 'EMMA',
        standard_name: '',
        standard_definition: null
      },
      queryStatus: {
        polling: {},
        result: {},
        requestError: ''
      },
      infos: infos
    }
  },
  infos: infos,
  methods: {
    validateForm () {
      var errors = []
      if (!this.form.parts) {
        errors.push('Provide parts !')
      }
      return errors
    }
  },
  watch: {
    'form.quantity_unit': function (val) {
      this.$set(this.form, 'part_quantity', this.quantityRanges[val].default)
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

.el-input-number.inline {
  margin-bottom: -9px;
  margin-left: 10px;
  width: 130px;
}

p.loadData {
  font-size: 0.8em;
  margin-bottom: 0;
  cursor: pointer;
}

.bands-range {
  .el-slider {
    display: inline-block;
    width: 150px;
    margin-bottom: -12px;
    margin-left: 10px;
  }
}
</style>
