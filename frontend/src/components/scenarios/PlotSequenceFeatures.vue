<template lang="pug">
.page
  //- el-alert(type='error' center :closable="false" title='Do not use !', show-icon)
  //-   .inline This app is currently under development
  h1  {{ infos.title }}
  img.icon.center-block(slot='title-img', :src='infos.icon')
  p.scenario-description Plot assembly features, highlight the important, discard the irrelevant, use the colors you love, etc.
  web-links(:emailSubject="'[CUBA] Feedback on web app: ' + infos.title",
            tweetMessage="Find common regions between different DNA sequences:",
            :tweetUrl="'http://cuba.genomefoundry.org/' + infos.path")

  .form

    h4.formlabel Upload sequence files
    filesuploader(v-model='form.files', text="Drop files (or click to select)",
                  help='Fasta, Genbank, or Snapgene files. No file too large please :)', :multiple='true')
    h4.formlabel Customize the plots
    el-tabs
      el-tab-pane(label='Plot parameters')
        p.inline Display:
          span.options
            el-radio(class='radio' v-model='form.display' label='linear') Linear
            el-radio(class='radio' v-model='form.display' label='circular') Circular
        p.inline.plot-size Width (inches)
          el-input-number.inline(v-model="form.plot_width", size="mini",
                                   :min='4', :max='16')
        p
          el-checkbox(v-model='form.plot_ruler') Plot ruler (location indices)
        span(v-if="form.display == 'linear'")

          p
            el-checkbox(v-model='form.plot_full_sequence') Plot the full sequence's span
          p.inline(v-if='!form.plot_full_sequence') Plot segment
            el-input-number.inline(v-model="form.plot_from_position", size="mini",
                                   :min='1', :max='form.plot_to_position-1')
            span &nbsp;to
            el-input-number.inline(v-model="form.plot_to_position", size="mini",
                                   :min='form.plot_from_position+1', :max='1000000')
          p
            el-checkbox(v-model='form.inline_labels') Allow inline feature labels.
          p
            el-checkbox(v-model='form.plot_nucleotides') Indicate sequence nucleotides
          p
            el-checkbox(v-model='form.plot_translation') Indicate amino acids
          p.inline(v-if='form.plot_translation') Plot segment
            el-input-number.inline(v-model="form.translation_start", size="mini",
                                   :min='1', :max='form.translation_end-3')
            span &nbsp;to
            el-input-number.inline(v-model="form.translation_end", size="mini",
                                   :min='form.translation_start+3', :max='1000000')
      el-tab-pane(label='Feature filters')
        p
          el-radio(v-model='form.keep_or_discard' label='keep') Only show these types
          el-radio(v-model='form.keep_or_discard' label='discard') Discard these types
        el-select(v-model='form.keep_or_discard_types' multiple placeholder="Select feature types")
          el-option(v-for='type, i in featureTypes', :label='type', :value='type', :key='type')
        p Features must contain
        el-input(v-model='form.keep_text',
                 placeholder='Enter comma-separated text fragment')
        p Features must NOT contain
        el-input(v-model='form.discard_text',
                 placeholder='Enter comma-separated text fragment')
      el-tab-pane(label='Feature styles')
        p.inline Default color:
          span.options
            el-color-picker.inline(v-model="form.default_color", size="small")
        p.inline Default line thickness:
          span.options
            el-radio(v-model="form.default_thickness", :label='0') 0
            el-radio(v-model="form.default_thickness", :label='1') 1
            el-radio(v-model="form.default_thickness", :label='2') 2
        p
          el-checkbox(v-model='form.default_display_label') Display labels by default

        el-card.style-definer(v-for='style, i in form.custom_styles', :key='i')
          p
            el-radio(v-model='style.keep_or_discard' label='keep') Features with
            el-radio(v-model='style.keep_or_discard' label='discard') Features without
          el-select(v-model='style.selector')
            el-option(label='Feature type' value='type')
            el-option(label='Feature text' value='text')
          el-select(v-if="style.selector === 'type'" v-model='style.feature_type')
            el-option(v-for='type, i in featureTypes', :label='type', :value='type', :key='type')
          el-input(v-if="style.selector === 'text'" v-model='style.feature_text',
                   placeholder='Enter an exact text fragment')
          p.inline Color:
            span.options
              el-color-picker.inline(v-model="style.color", size="small")
            span.options
              el-checkbox(v-model='style.display_label') Display label
          p.inline Line thickness:
            span.options
              el-radio(v-model="style.thickness", :label='0') 0
              el-radio(v-model="style.thickness", :label='1') 1
              el-radio(v-model="style.thickness", :label='2') 2

        el-button.center(icon='el-icon-plus' size='mini' @click='addCustomStyle()') Add a new style

    h4.formlabel RENDER THE PLOTS
    p.center
      el-radio(v-model='form.pdf_report', :label='false') as SVG images
      el-radio(v-model='form.pdf_report', :label='true') as a PDF

    backend-querier(:form='form', :backendUrl='infos.backendUrl',
                    :validateForm='validateForm', submitButtonText='Render the plots',
                    v-model='queryStatus')
    el-alert(v-show='queryStatus.requestError', :title="queryStatus.requestError",
       type="error", :closable="false")

  .results(v-if='!queryStatus.polling.inProgress && queryStatus.result')
    download-button(v-if='queryStatus.result.pdf_report',
                    :filedata='queryStatus.result.pdf_report')
    .figures-preview(v-if='queryStatus.result.figures_data')
      .figure-preview(v-for='fig in queryStatus.result.figures_data')
        h4 {{fig.filename}}
        img(:src='fig.img_data')


  powered-by(:softwareNames='infos.poweredby')
</template>

<script>
import learnmore from '../../components/widgets/LearnMore'
import filesuploader from '../../components/widgets/FilesUploader'

var infos = {
  title: 'Plot sequence features',
  navbarTitle: 'Plot sequence features',
  path: 'plot-sequence-features',
  description: '',
  backendUrl: 'start/plot_sequence_features',
  icon: require('assets/images/plot_sequence_features.svg'),
  poweredby: ['dnafeaturesviewer']
}

export default {
  data () {
    return {
      infos: infos,
      form: {
        files: [],
        display: 'linear',
        default_color: '#BABAEB',
        default_display_label: true,
        default_thickness: 1,
        plot_width: 12,
        plot_ruler: true,
        inline_labels: false,
        plot_full_sequence: true,
        plot_from_position: 1,
        plot_to_position: 10000,
        plot_nucleotides: false,
        plot_translation: false,
        translation_start: 0,
        translation_end: 10000,
        custom_styles: [],
        must_contain: '',
        must_not_contain: '',
        keep_or_discard: 'keep',
        keep_or_discard_types: [],
        pdf_report: false
      },
      default_style: {
        selector: 'type',
        keep_or_discard: 'keep',
        feature_type: 'CDS',
        feature_text: '',
        display_label: true,
        color: '#BABAEB',
        thickness: 1
      },
      queryStatus: {
        polling: {},
        result: {},
        requestError: ''
      },
      featureTypes: ['CDS', 'Promoter', 'Terminator', 'Source', 'Operon', 'misc_feature']
    }
  },
  components: {
    filesuploader,
    learnmore
  },
  infos: infos,
  methods: {
    validateForm () {
      var errors = []
      if (!this.form.files.length) {
        errors.push('Provide files !')
      }
      return errors
    },
    addCustomStyle () {
      this.form.custom_styles.push(Object.assign({}, this.default_style))
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

.form /deep/ {
  .el-color-picker__trigger {
    margin-bottom: -9px;
  }
  .el-checkbox__label {
    font-weight: normal;
    color: #2c3e50;
    font-size: 16px;
  }
}

.title-img {
  height:80px;
  margin-top: -20px;
  margin-bottom: 20px;
}



.el-checkbox.inline {
  margin-left: 15px;
}


.el-slider.inline {
  display: inline-block;
  margin-left: 20px;
  margin-right: 20px;
  margin-bottom: -12px;
}

.el-input-number.inline {
  margin-bottom: -9px;
  margin-left: 10px;
  width: 130px;
}

.plot-size {
  .el-input-number {
    width: 100px;
  }
}

.el-select {
  width: 100%
}

.result_image {
  max-width: 80%;
  margin: 0 10%;
}

.style-definer {
  // background-color: #fafafa;
  // padding: 10px;
  margin-bottom: 2em;
  width: 80%;
  margin-left: 10%
}

.figure-preview {
  margin-bottom: 5em;
  img {
    width:100%;
  }
}
</style>
