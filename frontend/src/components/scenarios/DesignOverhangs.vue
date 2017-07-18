<template lang="pug">
div
  h1  {{ infos.title }}
  img.icon.center-block(slot='title-img', :src='infos.icon')
  p.center Find sets of compatible overhangs for your assembly problem.
  learnmore Bla bla bla

  .form

    h4.formlabel What do you want to design ?

    el-select(v-model='form.goal', placeholder='Select')
      el-option(v-for='item in goal_options', :label='item.label',
                :value='item.value', :key='item.value')
    div(v-show='form.goal')
      div(v-show="form.goal === 'sequence_decomposition'")


        h4.formlabel Sequence


        filesuploader(v-model='form.sequence', text="Drop a single file (or click to select)",
                      help='Fasta or genbank. No file too large please :)', :multiple='false')
        p.file-uploaded(v-if='form.sequence') File selected: <b> {{form.sequence.name}} </b>


        h4.formlabel
          span Cutting mode
          helper(help="If the provided sequence has features with a label '#cut',\
                       these features can be used as cutting zones. Otherwise the\
                       sequence will be cut into similar-length fragments")

        p
          el-radio(class='radio' v-model='form.cutting_mode' label='equal') Cut equal fragments
          el-radio(class='radio' v-model='form.cutting_mode' label='features') Cut in featured zones

        p
          el-checkbox(v-model='form.extremities') Consider overhangs on extremities
          helper(help="If this is checked, the overhangs at the suggested cut\
                       positions will also be compatible with the 4 basepairs\
                       at each extremity of the sequence")

        p(v-if="form.cutting_mode === 'equal'")
          span Number of fragments:
          el-input-number.inline(v-model="form.n_fragments", size="small",
                                 :min=2, :max=60)





      h4.formlabel Overhangs parameters

      p.inline(v-if="form.goal === 'overhangs_set'")
        span Number of overhangs:
        el-input-number.inline(v-show='!form.auto_overhangs'
                               v-model="form.n_overhangs", size="small",
                               :min=1, :max=50)
        el-checkbox.inline(v-model='form.auto_overhangs') as many as possible

      //- p.inline
      //-   span Overhang size, in nucleotides:
      //-   el-input-number.inline(v-model="form.overhang_size", size="small",
      //-                          :min=2, :max=6)

      p.inline Differences between overhangs:
        el-input-number.inline(v-model="form.overhangs_differences", size="small",
                               :min=1, :max=4)
        helper(help='Number of basepairs by which every overhang pair should\
                     differ (this includes reverse complements).')

      p.inline
        span Overhangs GC content
        el-slider.inline(range show-stops v-model='form.gc_content',
                         :max=100, :min=0, :step=25, :style="{width: '100px'}")
        span {{ form.gc_content[0] }}% to {{ form.gc_content[1] }}%:

      h4 Forbidden overhangs

      el-select(v-model='form.forbidden_overhangs',
                placeholder='Enter any forbidden overhangs e.g. ATGC, GCTA, ...',
                filterable, :multiple='true')
        el-option(v-for='overhang in overhangs', :label='overhang', :value='overhang', :key='overhang')

      div(v-show="form.goal === 'overhangs_set'")
        h4 Mandatory overhangs
        el-select(v-model='form.mandatory_overhangs',
                  placeholder='Enter any mandatory overhangs e.g. ATGC, GCTA, ...',
                  filterable, :multiple='true')
          el-option(v-for='overhang in overhangs', :label='overhang', :value='overhang', :key='overhang')




    backend-querier(:form='form', :backendUrl='infos.backendUrl',
                    :validateForm='validateForm', submitButtonText='Evaluate',
                    v-show='form.goal', v-model='queryStatus')
    p(v-if='queryStatus.polling.inProgress && queryStatus.polling.data.n_overhangs').center Attempting to find {{queryStatus.polling.data.n_overhangs}} overhangs
    el-alert(v-show='queryStatus.requestError', :title="queryStatus.requestError",
       type="error", :closable="false")

  .results(v-if='!queryStatus.polling.inProgress && queryStatus.polling.data')
    p(v-if='!queryStatus.result.success') No solution found 😢. Maybe your parameters are too restrictive ?
    div(v-if='selected_overhangs')
      p We found a collection of {{ selected_overhangs.length }} overhangs:
      p.selected-overhangs
        span(v-for='overhang in selected_overhangs', :key='overhang') {{overhang}}
    download-button(v-if='queryStatus.result.zip_file',
                    :filedata='queryStatus.result.zip_file')

  powered-by(:softwareNames='infos.poweredby')
</template>

<script>
import learnmore from '../../components/widgets/LearnMore'
import filesuploader from '../../components/widgets/FilesUploader'
import tools from '../../tools'

var infos = {
  title: 'Design Golden Gate Overhangs',
  navbarTitle: 'Design Golden Gate Overhangs',
  path: 'design-overhangs',
  description: '',
  backendUrl: 'start/design_overhangs',
  icon: require('assets/images/overhangs.svg'),
  poweredby: ['goldenhinges', 'dnachisel']
}

export default {
  data: function () {
    return {
      overhangs: tools.cartesian(['A', 'T', 'G', 'C'],
                                 ['A', 'T', 'G', 'C'],
                                 ['A', 'T', 'G', 'C'],
                                 ['A', 'T', 'G', 'C']).map((l) => l.join('')),
      infos: infos,
      form: {
        goal: null,
        overhangs_differences: 1,
        gc_content: [25, 75],
        sequence: null,
        left_flank_sequence: '',
        right_flank_sequence: '',
        mandatory_overhangs: [],
        forbidden_overhangs: [],
        cutting_mode: 'equal',
        n_overhangs: 10,
        n_fragments: 1,
        auto_overhangs: true,
        extremities: true
      },
      queryStatus: {
        polling: {},
        result: {},
        requestError: ''
      },
      goal_options: [
        {
          label: 'A collection of compatible overhangs',
          value: 'overhangs_set'
        },
        {
          label: 'A sequence decomposition, with compatible overhangs',
          value: 'sequence_decomposition'
        }
      ]
    }
  },
  components: {
    filesuploader,
    learnmore
  },
  infos: infos,
  methods: {
    handleSuccess: function (evt) {
      console.log(evt)
    },
    validateForm: function () {
      var errors = []
      if (!this.form.goal) {
        errors.push('Provide a goal !')
      } else if ((this.form.goal === 'sequence_decomposition') && (!this.form.sequence)) {
        errors.push('Provide a sequence !')
      }
      return errors
    }
  },
  computed: {
    selected_overhangs: function () {
      if (!this.queryStatus.result) {
        return null
      } else {
        return this.queryStatus.result.overhangs
      }
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
  font-size: 16px !important;
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

.file-uploaded {
  font-size: 14px;
  margin-top: -10px;
}
</style>
