<template lang="pug">
.page
  h1 {{infos.title}}
  web-links(:emailSubject="'[CUBA] Feedback on web app: ' + infos.title",
            tweetMessage="Get enzyme suggestions for your restriction digests:",
            :tweetUrl="'https://cuba.genomefoundry.org/' + infos.path")
  img.icon.center-block(slot='title-img', :src='infos.icon')
  p.scenario-description.
    Provide a source plate map and an assembly plan, get a robotic picklist
    spreadsheet for Tecan EVO or Labcyte Echo.

  el-alert(title="This app is a stub. Not all features may work properly" type="warning" show-icon)

  .form
    h4.formlabel Assembly Plan
    collapsible(title='Examples')
      file-example(filename='example_picklist.xls',
                   @input='function (e) {form.picklist = e}',
                   fileHref='/static/file_examples/create_assembly_picklists/example_picklist.xls',
                   imgSrc='/static/file_examples/generic_logos/spreadsheet.svg')
        p.
          Picklist to build 7 assemblies using common parts of the EMMA standard

    files-uploader(v-model='form.picklist', tip='Excel file or CSV', :multiple='false')
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

    hr
    h4.formlabel Mix
    el-form(label-width='200px')
      el-form-item(label='Part quantity unit')
        el-select(v-model='form.quantity_unit')
          el-option(label='femto-mole', value='fmol')
          el-option(label='nano-gram', value='ng')
          el-option(label='nano-molar', value='nM')

      el-form-item(:label="`Part quantity (${form.quantity_unit})`")
        el-input-number(v-model='form.part_quantity',
                        :min='quantityRanges[form.quantity_unit].min',
                        :max='quantityRanges[form.quantity_unit].max'
                        :step='quantityRanges[form.quantity_unit].step')
      el-form-item(label='Buffer volume (µL)')
        el-input-number(v-model='form.buffer_volume', :min='0', :max='5000', :step='0.05')
      el-form-item(label='Total volume (µL)')
        el-input-number(v-model='form.total_volume', :min='0', :max='2000', :step='0.05')


      //- el-checkbox(v-if="form.quantity_unit !== 'ng'", v-model='form.specify_ratio') 

      p(v-if="form.quantity_unit !== 'ng'")
        | Specify part:backbone molar ratio. The plasmid backbone name must match the 
        | name in the source plate (case-sensitive).
        | You may specify multiple backbones, separated by commas.
      el-form-item(v-if="form.quantity_unit !== 'ng'", label='Part : backbone ratio')
        el-input-number(v-model='form.part_backbone_ratio', :min='0.01', :max='100', :step='0.5')
      el-form-item(v-if="form.quantity_unit !== 'ng'", label='Backbone name')
        el-input(v-model='form.backbone_name')

    .parts-infos(v-if="form.quantity_unit !== 'ng'")
      hr
      h4.formlabel Part information
      :markdown-it
        Since a molar quantity is desired, we need to know the molecular weight of the DNA.
        Provide either
        - A table of the part sizes, backbone included (see example file) or
        - A series of Genbank/Fasta sequences of the parts.

      collapsible(title='Examples')
        file-example(filename='emma_parts.zip',
                     @input='function (e) {form.parts_infos.push(e)}',
                     fileHref='/static/file_examples/create_assembly_picklists/emma_parts.zip',
                     imgSrc='/static/file_examples/generic_logos/sequences_records.png')
          p.
            ZIP archive containing the sequences of standard EMMA parts.
        file-example(filename='emma_parts_sizes.csv',
                     @input='function (e) {form.parts_infos.push(e)}',
                     fileHref='/static/file_examples/create_assembly_picklists/emma_parts_sizes.csv',
                     imgSrc='/static/file_examples/generic_logos/spreadsheet.svg')
          p.
            Spreadsheet indicating the total size of EMMA parts (including vector backbone).

      files-uploader(v-model='form.parts_infos', tip='Genbank/Fasta/Snapgene sequences, or spreadsheet')
      el-checkbox(v-model='form.use_file_names_as_ids') Use file names (not record IDs) to identify parts.
    hr
    h4.formlabel Source plate
    :markdown-it
      The volume should be specified in µL and the concentration in ng/µL (see example file).

    collapsible(title='Examples')
      file-example(filename='example_echo_plate.xlsx',
                   @input='function (e) {form.source_plate =e}',
                   fileHref='/static/file_examples/create_assembly_picklists/example_echo_plate.xlsx',
                   imgSrc='/static/file_examples/generic_logos/spreadsheet.svg')
        p.
          Map of a plate with common genetic parts from the EMMA standard.

    files-uploader(v-model='form.source_plate', tip='Excel file with "content", "volume" and "concentration" sheets', :multiple='false')
    hr
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

      collapsible(title='Examples')
        file-example(filename='example_destination_plate.xlsx',
                     @input='function (e) {form.destination_plate = e}',
                     fileHref='/static/file_examples/create_assembly_picklists/example_destination_plate.xlsx',
                     imgSrc='/static/file_examples/generic_logos/spreadsheet.svg')
          p.
            Map of a plate with common genetic parts from the EMMA standard
      files-uploader(v-model='form.destination_plate', :multiple='false')

    backend-querier(:form='form', :backendUrl='infos.backendUrl',
                    :validateForm='validateForm', submitButtonText='Generate a picklist',
                    v-model='queryStatus')
    progress-bars(:bars='queryStatus.polling.data.bars', :order="['record', 'primer']"
                  v-if='queryStatus.polling.inProgress && queryStatus.polling.data')
    el-alert(v-if='queryStatus.requestError', :title="queryStatus.requestError",
       type="error", :closable="false")

    .results(v-if='!queryStatus.polling.inProgress')
      p(v-if='queryStatus.result.message') {{ queryStatus.result.message }}
      div(v-if='queryStatus.result.missing_parts')
        ul
          li(v-for='didYouMean, part in queryStatus.result.missing_parts')
            span #[b {{part}}]. <br/>
            span(v-if='didYouMean.featured_in') Featured in {{didYouMean.featured_in.join(', ')}} <br/>
            span Did you mean: <br/>
            ul(v-if='didYouMean.featured_in')
              li(v-for='part in didYouMean.did_you_mean') {{part[0]}}, {{part[1]}}
            ul(v-else)
              li(v-for='part in didYouMean') {{part}}
      div(v-if='queryStatus.result.picklist_data && queryStatus.result.picklist_data.duplicates'
          style='color: red')
        ul
          li(v-for='duplicates, part_name in queryStatus.result.picklist_data.duplicates').
            Duplicates wells for #[b {{part_name}}]: {{duplicates.wells}}
            ({{duplicates.selected}} was selected)
      download-button(v-if='queryStatus.result.file',
        text='Download Picklist',
        :filedata='queryStatus.result.file')
  powered-by(:softwareNames='infos.poweredby')
</template>

<script>
var infos = {
  title: "Create Assembly Picklists",
  navbarTitle: "Create Assembly Picklists",
  path: "create_assembly_picklists",
  description: "",
  backendUrl: "start/create_assembly_picklists",
  icon: require("../../assets/images/create_assembly_picklists.svg"),
  poweredby: ["plateo"],
};

export default {
  data() {
    return {
      form: {
        picklist: null,
        source_plate: null,
        destination_type: "new",
        destination_size: 96,
        destination_plate: null,
        fill_by: "column",
        quantity_unit: "fmol",
        part_quantity: 1.3,
        buffer_volume: 0.3,
        total_volume: 1,
        part_backbone_ratio: 1,
        backbone_name: "hc_amp_ccdb",
        parts_infos: [],
        dispenser_machine: "labcyte_echo",
        dispenser_max_volume: 0.5,
        dispenser_min_volume: 5,
        dispenser_resolution: 2.5,
        dispenser_dead_volume: 8,
        use_file_names_as_ids: true,
      },
      quantityRanges: {
        fmol: {
          min: 0.1,
          max: 100,
          step: 0.1,
          default: 1.3,
        },
        nM: {
          min: 0.1,
          max: 10,
          step: 0.1,
          default: 1.3,
        },
        ng: {
          min: 1,
          max: 1000,
          step: 0.1,
          default: 10,
        },
      },
      queryStatus: {
        polling: {},
        result: {},
        requestError: "",
      },
      infos: infos,
    };
  },
  infos: infos,
  methods: {
    validateForm() {
      var errors = [];
      if (!this.form.source_plate) {
        errors.push("Provide a source plate !");
      }
      return errors;
    },
  },
  watch: {
    "form.quantity_unit": function (val) {
      this.$set(this.form, "part_quantity", this.quantityRanges[val].default);
    },
  },
};
</script>

<style lang='scss' scoped>
h4.formlabel {
  text-align: center;
  text-transform: uppercase;
  margin-top: 40px;
}

.form {
  margin: 50px auto;
  max-width: 500px;
}

.title-img {
  height: 80px;
  margin-top: -20px;
  margin-bottom: 20px;
}

.el-checkbox {
  font-weight: normal;
}

.el-select {
  width: 100%;
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
