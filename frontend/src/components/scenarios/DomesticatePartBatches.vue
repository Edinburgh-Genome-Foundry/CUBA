<template lang="pug">
.page
  h1 {{infos.title}}
  web-links(:emailSubject="'[CUBA] Feedback on web app: ' + infos.title",
            tweetMessage="Domesticate parts for EMMA and other standards:",
            :tweetUrl="'https://cuba.genomefoundry.org/' + infos.path")
  img.icon.center-block(slot='title-img', :src='infos.icon')
  p.scenario-description.
    Provide a part sequences to domesticate for EMMA assembly or another standard.

  .form
    h4.formlabel Standard

    //- p
    //-   el-select(v-model='form.standard')
    //-     el-option(value='EMMA', label='EMMA')
    //-     el-option(value='custom', label='Custom')
    el-form
      collapsible(title='Examples')
        file-example(filename='EMMA_standard.csv',
                      @input="function (e) { \
                        form.standard_definition = e; \
                        form.standard_name = 'EMMA'; \
                      }",
                      fileHref='/static/file_examples/domesticate_part_batches/EMMA.csv',
                      imgSrc='/static/file_examples/generic_logos/spreadsheet.svg')
          p.
            A table representing the EMMA standard. This standard requires
            2-nucleotide additions in some slots for the parts to follow
            in frame: the two nucleotides and the 4-nucleotide overhang
            create a 6-nucleotide scar (2 amino acids in the translated
            sequence).
      el-form-item(label='Name')
        el-input(v-model='form.standard_name', placeholder='Name of the standard')
      el-form-item(label='Definition')
        files-uploader(v-model='form.standard_definition', tip='csv or excel', :multiple='false')
    hr
    h4.formlabel Parts
    collapsible(title='Examples')
      file-example(filename='EMMA_undomesticated_parts.zip',
                    @input="function (e) { \
                      form.parts = [e]; \
                    }",
                    fileHref='/static/file_examples/domesticate_part_batches/EMMA_undomesticated_parts.zip',
                    imgSrc='/static/file_examples/generic_logos/linear_part_records.svg')
        p.
          Linear records of parts to be domesticated for the EMMA standard
    files-uploader(v-model='form.parts', tip='Genbank/Fasta/Snapgene sequences, or spreadsheet')
    el-select(v-model='form.part_id_source')
      el-option(value='file_name' label='Use file names (without extensions) as part IDs')
      el-option(value='record_id' label='Use Genbank/Fasta record ID/name as part IDs')
    hr
    el-checkbox(v-model='form.allow_edits') Allow sequence edits
    
    backend-querier(:form='form', :backendUrl='infos.backendUrl',
                    :validateForm='validateForm', submitButtonText='Domesticate',
                    v-model='queryStatus')
    progress-bars(:bars='queryStatus.polling.data.bars', :order="['record', 'primer']"
                  v-if='queryStatus.polling.inProgress && queryStatus.polling.data')
    el-alert(v-if='queryStatus.requestError', :title="queryStatus.requestError",
       type="error", :closable="false")

    .results(v-if='!queryStatus.polling.inProgress && queryStatus.polling.data')
      div(v-if='queryStatus.result.nfails > 0')
        p.
          There were {{queryStatus.result.nfails}} errors, see report for more.
        p(v-if='!form.allow_edits').
          Maybe you should tick "Allow sequence edits" so sequences can be mutated to fix errors.
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
        part_id_source: 'record_id',
        standard: 'custom',
        standard_name: '',
        allow_edits: false,
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
