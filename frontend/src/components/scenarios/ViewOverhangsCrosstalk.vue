<template lang="pug">
.page
  h1  {{ infos.title }}
  img.icon.center-block(slot='title-img', :src='infos.icon')
  p.scenario-description Find sets of compatible overhangs for your assembly problem.
  web-links(:emailSubject="'[CUBA] Feedback on web app: ' + infos.title",
            tweetMessage="Browse DNA overhangs annealing and cross-talk data:",
            :tweetUrl="'https://cuba.genomefoundry.org/' + infos.path")
  //- learnmore Bla bla bla

  .form

    h4.formlabel Your overhangs
    .example(@click="form.overhangs = 'ATTG, GCTA, ACAT, CCAG'"
             style='cursor: pointer; color: grey; margin: 0.5em;') Example
    el-input(type='textarea',
             :rows='4',
             v-model='form.overhangs',
             placeholder='Enter comma-separated overhangs, e.g. ATTG, TTAC, ...')
    p
      .dataset-parameter Dataset:
      .dataset-parameter
        el-select(v-model='form.temperature')
          el-option(value='25C', label='25C')
          el-option(value='37C', label='37C')
      .dataset-parameter
        el-select(v-model='form.incubation')
          el-option(value='01h', label='1h')
          el-option(value='18h', label='18h')
        


    backend-querier(:form='form', :backendUrl='infos.backendUrl',
                    :validateForm='validateForm', submitButtonText='Plot crosstalk',
                    v-model='queryStatus')
    el-alert(v-show='queryStatus.requestError', :title="queryStatus.requestError",
       type="error", :closable="false")

  .results(v-if='!queryStatus.polling.inProgress && queryStatus.result.figure_data')
    center
      img.result_image(:src='queryStatus.result.figure_data')

    p(style="margin-bottom: 5em").
      In the this plot, if you see anything else than the square pairs around
      the diagonal, it means there is cross-talking between your overhangs
      (so risk of misannealing). If one of these diagmonal square pairs appears
      lighter than the others, it means that the corresponding overhang has
      weak self-annealing (risk of having no assembly).
  
  :markdown-it
    **About this app**

    This app downloads experimental data from
    [this BioRxiv paper](https://www.biorxiv.org/content/early/2018/05/15/322297),
    *Optimization of Golden Gate assembly through application of ligation
    sequence-dependent fidelity and bias profiling*, Potapov *et al.*, 2018.

    It reads the spreadsheet provided in supplementary information, keeps only
    the data relevant to your overhangs, and plots it as a heatmap, similar to
    the figures shown in the paper.

    




  powered-by(:softwareNames='infos.poweredby')
</template>

<script>
import learnmore from '../../components/widgets/LearnMore'

var infos = {
  title: 'View Overhangs Crosstalk',
  navbarTitle: 'View Overhangs Crosstalk',
  path: 'view-overhangs-crosstalk',
  description: '',
  backendUrl: 'start/view_overhangs_crosstalk',
  icon: require('../../assets/images/view_overhangs_crosstalk.svg'),
  poweredby: ['tatapov']
}

export default {
  data () {
    return {
      infos: infos,
      form: {
        overhangs: '',
        temperature: '37C',
        incubation: '01h'
      },
      queryStatus: {
        polling: {},
        result: {},
        requestError: ''
      }
    }
  },
  components: {
    learnmore
  },
  infos: infos,
  methods: {
    validateForm () {
      var errors = []
      if (!this.form.overhangs.length) {
        errors.push('Provide files !')
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


.result_image {
  max-width: 80%;
  margin: 0 10%;
}

.dataset-parameter {
  display: inline-block;
  width: 20%;
  margin-left: 2em;
}
</style>
