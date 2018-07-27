<template lang="pug">
.page
  h1 {{infos.title}}
  web-links(:emailSubject="'[CUBA] Feedback on web app: ' + infos.title",
            :tweetMessage="infos.tweetMessage",
            :tweetUrl="'http://cuba.genomefoundry.org/' + infos.path")
  img.icon.center-block(slot='title-img', :src='infos.icon' title='NOT a proto-nazi symbol !')
  p.scenario-description.
    Provide a source plate map and a list of constructs, get a robotic picklist (=spreadsheet)
    for Tecan EVO or Labcyte Echo.

  el-alert(title="This is a very young app, no warranties." type="warning" show-icon)

  .form
    h4.formlabel Source Plate
    collapsible(title='Examples')
      file-example(filename='example_source_plate.xlsx',
                   @input='function (e) {form.source_plate = e}',
                   fileHref='/static/file_examples/rearray_plates/example_rearray_source_plate.xlsx',
                   imgSrc='/static/file_examples/generic_logos/spreadsheet.svg')
        p.
          Picklist to build 7 assemblies using common parts of the EMMA standard
    files-uploader(v-model='form.source_plate', tip='Excel file or CSV', :multiple='false')
    el-form(label-width='150px')
      el-form-item(label='Sample by')
        el-select(v-model='form.samples_source_by')
          el-option(value='row' label='Row')
          el-option(value='column' label='Column')
    hr
    h4.formlabel Destination Plate
    el-form(label-width='150px')
      el-form-item(label='Plate')
        el-select(v-model='form.destination_type')
          el-option(value='new', label='New Plate')
          el-option(value='file', label='Existing Plate')
      el-form-item(label='Plate size' v-if="form.destination_type === 'new'")
        el-select(v-model='form.destination_size')
          el-option(:value='96', label='96 wells')
          el-option(:value='384', label='384 wells')
          el-option(:value='1536', label='1536 wells')
      el-form-item(label='Destination map' v-if="form.destination_type === 'file'")
        files-uploader(v-model='form.destination_plate', :multiple='false')
        collapsible(title='Examples')
          file-example(filename='example_destination_plate.xlsx',
                         @input='function (e) {form.destination_plate = e}',
                         fileHref='/static/file_examples/rearray_plates/example_destination_plate.xlsx',
                         imgSrc='/static/file_examples/generic_logos/spreadsheet.svg')
              p.
                Map of a plate with common genetic parts from the EMMA standard

      el-form-item(label='Fill by')
        el-select(v-model='form.fill_by')
          el-option(value='row' label='Row')
          el-option(value='column' label='Column')
    hr
    h4.formlabel Rearraying
    el-form(label-width='250px')
      el-form-item(label='Volume per well (µL)')
        el-select(v-model='form.total_volume')
      el-form-item(label='Normalize concentrations')
        el-checkbox(v-model='form.normalize')
      el-form-item(v-if='form.normalize' label='Concentration unit')
        el-select(v-model='form.concentration_unit')
          el-option(label='ng/µL', value='ngul')
          el-option(label='fM', value='fM')
      el-form-item(v-if='form.normalize', :label="'Concentration (' + concentrationUnit +')'")
        el-input-number(v-model='form.concentration')
    hr
    h4.formlabel Dispenser
    el-form(label-width='250px')
      el-form-item(label='Machine')
        el-select(v-model='form.dispenser_machine')
          el-option(label='Tecan EVO', value='tecan_evo')
          el-option(label='Labcyte ECHO', value='labcyte_echo')
      el-form-item(label='Min volume (nL)')
        el-input-number(v-model='form.dispenser_min_volume', :min='0', :max='1000', :step='0.5')
      el-form-item(label='Max one-time volume (µL)')
        el-input-number(v-model='form.dispenser_max_volume', :min='0', :max='2000', :step='0.5')
      el-form-item(label='Volume resolution (nL)')
        el-input-number(v-model='form.dispenser_resolution', :min='0', :max='1000', :step='0.5')
      el-form-item(label='Dead volume (µL)')
        el-input-number(v-model='form.dispenser_dead_volume', :min='0', :max='1000', :step='0.1')

    .parts-infos(v-if="form.normalize && (form.concentration_unit === 'fM')")
      hr
      h4.formlabel Parts infos
      :markdown-it
        Since a molar quantity is desired, we need to know the molecular weigth of the DNA.
        Provide either
        - A table of the parts sizes, backbone included (see example file) or
        - A series of Genbank/Fasta sequences of the parts.

      collapsible(title='Examples')
        file-example(filename='emma_parts.zip',
                     @input='function (e) {form.parts_infos.push(e)}',
                     fileHref='/static/file_examples/rearray_plates/emma_parts.zip',
                     imgSrc='/static/file_examples/generic_logos/sequences_records.png')
          p.
            ZIP archive containing the sequences of standard EMMA parts.
        file-example(filename='emma_parts_sizes.csv',
                     @input='function (e) {form.parts_infos.push(e)}',
                     fileHref='/static/file_examples/rearray_plates/emma_parts_sizes.csv',
                     imgSrc='/static/file_examples/generic_logos/spreadsheet.svg')
          p.
            Spreadsheet indicating the total size of EMMA parts (including vector backbone).

      files-uploader(v-model='form.parts_infos', tip='Genbank/Fasta/Snapgene sequences, or spreadsheet')

    backend-querier(:form='form', :backendUrl='infos.backendUrl',
                    :validateForm='validateForm', submitButtonText='Generate a picklist',
                    v-model='queryStatus')
    progress-bars(:bars='queryStatus.polling.data.bars', :order="['record', 'primer']"
                  v-if='queryStatus.polling.inProgress && queryStatus.polling.data')
    el-alert(v-if='queryStatus.requestError', :title="queryStatus.requestError",
       type="error", :closable="false")

    .results(v-if='!queryStatus.polling.inProgress')
      download-button(v-if='queryStatus.result.file',
        text='Download Picklist',
        :filedata='queryStatus.result.file')
  powered-by(:softwareNames='infos.poweredby')
</template>

<script>

var infos = {
  title: 'Rearray Plates',
  navbarTitle: 'Rearray Plates',
  path: 'rearray_plates',
  description: '',
  tweetMessage: 'Create picklists to rearray plates (TECAN, EVO)',
  backendUrl: 'start/rearray_plates',
  icon: require('../../assets/images/rearray_plates2.svg'),
  poweredby: ['plateo']
}

export default {
  data () {
    return {
      form: {
        picklist: null,
        source_plate: null,
        destination_type: 'new',
        destination_size: 96,
        destination_plate: null,
        sample_source_by: 'column',
        fill_by: 'column',
        normalize: true,
        concentration_unit: 'ngul',
        concentration: 100,
        total_volume: 50,
        parts_infos: [],
        dispenser_machine: 'labcyte_echo',
        dispenser_max_volume: 0.5,
        dispenser_min_volume: 5,
        dispenser_resolution: 2.5,
        dispenser_dead_volume: 8
      },
      concentrationRanges: {
        'fM': {
          min: 0.1,
          max: 100,
          step: 0.1,
          default: 1.3
        },
        'ngul': {
          min: 1,
          max: 500,
          step: 0.1,
          default: 100
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
  computed: {
    concentrationUnit () {
      return (this.form.concentration_unit === 'ngul') ? 'ng/µL' : 'fM'
    }
  },
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
