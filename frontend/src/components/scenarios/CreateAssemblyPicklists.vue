<template lang="pug">
.page
  h1 {{infos.title}}
  web-links(:emailSubject="'[CUBA] Feedback on web app: ' + infos.title",
            tweetMessage="Get enzyme suggestions for your restriction digests:",
            :tweetUrl="'http://cuba.genomefoundry.org/' + infos.path")
  img.icon.center-block(slot='title-img', :src='infos.icon' title='NOT a proto-nazi symbol !')
  p.scenario-description.
    Provide a source plate map and a list of constructs, get a robotic picklist (=spreadsheet)
    for Tecan EVO or Labcyte Echo.

  .form
    h4.formlabel Parameters
    el-form(label-width='200px')
      el-form-item(label='Part quantity unit')
        el-select(v-model='form.quantity_unit')
          el-option(label='femtomoles', value='fM')
          el-option(label='nanograms', value='ng')

      el-form-item(:label="`Part quantity (${form.quantity_unit})`")
        el-input-number(v-model='form.part_quantity',
                        :min='quantityRanges[form.quantity_unit].min',
                        :max='quantityRanges[form.quantity_unit].max'
                        :step='quantityRanges[form.quantity_unit].step')
      el-form-item(label='Reagents volume (µL)')
        el-input-number(v-model='form.reagents_volume', :min='0', :max='5000', :step='0.05')
      el-form-item(label='Total volume (µL)')
        el-input-number(v-model='form.total_volume', :min='0', :max='2000', :step='0.05')
      el-form-item(label='Dispenser')
        el-select(v-model='form.dispenser')
          el-option(label='Tecan EVO', value='tecan_evo')
          el-option(label='Labcyte ECHO', value='labcyte_echo')
    .parts-infos(v-if="form.quantity_unit === 'fM'")
      h4.formlabel Parts infos
      p.
        Since a molar quantity is desired, we need to know the molecular weigth of the DNA.
        Provide either a table of the parts sizes, backbone included (see example file) or
        a series of Genbank/Fasta sequences of the parts.
      collapsible(title='Examples')
        file-example(filename='parts_records.zip',
                     @input='function (e) {form.source_plate.push(e)}',
                     fileHref='/static/file_examples/create_assembly_picklists/source_plate.xlsx',
                     imgSrc='/static/file_examples/generic_logos/spreadsheet.svg')
          p.
            Map of a plate with common genetic parts from the EMMA standard.
        file-example(filename='parts_sizes.xls',
                     @input='function (e) {form.source_plate.push(e)}',
                     fileHref='/static/file_examples/create_assembly_picklists/source_plate.xlsx',
                     imgSrc='/static/file_examples/generic_logos/spreadsheet.svg')
          p.
            Map of a plate with common genetic parts from the EMMA standard.

      files-uploader(v-model='form.parts_infos', tip='Genbank/Fasta/Snapgene sequences, or spreadsheet')
    h4.formlabel Source Plate
    collapsible(title='Examples')
      file-example(filename='source_plate.xlsx',
                   @input='function (e) {form.source_plate.push(e)}',
                   fileHref='/static/file_examples/create_assembly_picklists/source_plate.xlsx',
                   imgSrc='/static/file_examples/generic_logos/spreadsheet.svg')
        p.
          Map of a plate with common genetic parts from the EMMA standard.

    files-uploader(v-model='form.source_plate', tip='Fasta or Genbank files', :multiple='false')

    h4.formlabel Destination Plate
    el-form(label-width='150px')
      el-form-item(label='Plate')
        el-select(v-model='form.destination_type')
          el-option(value='new', label='New Plate')
          el-option(value='file', label='Existing Plate')
      el-form-item(label='Plate size')
        el-select(v-model='form.destination_size')
          el-option(:value='96', label='96 wells')
          el-option(:value='384', label='384 wells')
          el-option(:value='1536', label='1536 wells')
      el-form-item(label='Fill by')
        el-select(v-model='form.fill_by')
          el-option(value='row' label='Row')
          el-option(value='column' label='Column')
    .destination-file(v-if="form.destination_type === 'file'")
      files-uploader(v-model='form.destination_plate', :multiple='false')
      collapsible(title='Examples')
        file-example(filename='source_plate.xlsx',
                     @input='function (e) {form.source_plate.push(e)}',
                     fileHref='/static/file_examples/create_assembly_picklists/source_plate.xlsx',
                     imgSrc='/static/file_examples/generic_logos/spreadsheet.svg')
          p.
            Map of a plate with common genetic parts from the EMMA standard

    backend-querier(:form='form', :backendUrl='infos.backendUrl',
                    :validateForm='validateForm', submitButtonText='Select primers',
                    v-model='queryStatus')
    progress-bars(:bars='queryStatus.polling.data.bars', :order="['record', 'primer']"
                  v-if='queryStatus.polling.inProgress && queryStatus.polling.data')
    el-alert(v-if='queryStatus.requestError', :title="queryStatus.requestError",
       type="error", :closable="false")

    .results(v-if='!queryStatus.polling.inProgress')
      p.results-summary(v-if='queryStatus.result.summary',
                        v-html="queryStatus.result.summary")
      download-button(v-if='queryStatus.result.zip_file',
                      :filedata='queryStatus.result.zip_file')
  powered-by(:softwareNames='infos.poweredby')
</template>

<script>

var infos = {
  title: 'Create assembly picklists',
  navbarTitle: 'Create assembly picklists',
  path: 'create_assembly_picklists',
  description: '',
  backendUrl: 'start/create_assembly_picklists',
  icon: require('../../assets/images/create_assembly_picklists.svg'),
  poweredby: ['plateo']
}

export default {
  data () {
    return {
      form: {
        source_plate: null,
        destination_type: 'new',
        destination_size: 96,
        fill_by: 'column',
        destination_plate: null,
        quantity_unit: 'fM',
        part_quantity: 1.5,
        reagents_volume: 0.2,
        total_volume: 1,
        parts_infos: [],
        dispenser: 'labcyte_echo'
      },
      quantityRanges: {
        'fM': {
          min: 0.1,
          max: 100,
          step: 0.1,
          default: 1.5
        },
        'ng': {
          min: 1,
          max: 1000,
          step: 0.1,
          default: 10
        }
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
      if (!this.form.source_plate) {
        errors.push('Provide a source plate !')
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
